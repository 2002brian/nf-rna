"""Hardware-aware local resource policy for Nextflow scheduling.

This module decides only *how much of this machine* a run may schedule on.  It
never changes a tool argument that affects results: the aggregate ceiling
bounds Nextflow's local executor, and the single per-process adjustment
(Salmon quantification memory) is a scheduling request, not a Salmon option.

Detection is best effort.  Any failure degrades to the historical fixed
8 CPU / 12 GiB default instead of blocking a run.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

GIB = 1024**3
AUTO = "auto"

# The largest enabled local process contract (LARGE) needs 8 CPUs / 12 GiB.
# Auto mode never selects less than this when the machine itself provides it,
# so a small host behaves exactly as the former fixed default did.
FLOOR_CPUS = 8
FLOOR_MEMORY_GIB = 12
FALLBACK_CPUS = 8
FALLBACK_MEMORY_GIB = 12

# Conservative OS reserve: at least 1 CPU / 4 GiB, or 1/8 of CPUs / 15% of memory.
CPU_RESERVE_FRACTION = 0.125
CPU_RESERVE_MINIMUM = 1
MEMORY_RESERVE_FRACTION = 0.15
MEMORY_RESERVE_MINIMUM_GIB = 4

# cgroup v1 reports "no limit" as a huge page-aligned number.
_UNLIMITED_BYTES = 1 << 60

# nf-core/rnaseq 3.26.0 SALMON_QUANT is label process_medium (6 CPUs, 36 GB).
# Its peak RSS is dominated by the in-memory index, so the scheduling request
# is sized from the index on disk.  CPUs (Salmon --threads) are never changed.
SALMON_QUANT_SELECTOR = ".*:SALMON_QUANT"
SALMON_QUANT_NFCORE_MEMORY_GIB = 36
SALMON_INDEX_MEMORY_FACTOR = 1.15
SALMON_INDEX_MEMORY_OVERHEAD_GIB = 2


@dataclass(frozen=True)
class DetectedResources:
    """What this process may use, from the kernel's point of view.

    ``usable_*`` already respect CPU affinity and cgroup limits; the raw
    observations are kept in ``details`` for provenance.
    """

    usable_cpus: int | None = None
    usable_memory_bytes: int | None = None
    available_memory_bytes: int | None = None
    details: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        return {
            "usable_cpus": self.usable_cpus,
            "usable_memory_gib": _gib(self.usable_memory_bytes),
            "available_memory_gib": _gib(self.available_memory_bytes),
            **self.details,
        }


@dataclass(frozen=True)
class ProcessLimits:
    """Linux process constraints beyond the machine size; each is None when absent."""

    affinity_cpus: int | None = None
    cgroup_cpu_limit: float | None = None
    cgroup_memory_limit_bytes: int | None = None

    def apply_cpus(self, cpus: int | None) -> int | None:
        candidates = [value for value in (cpus, self.affinity_cpus) if value]
        if self.cgroup_cpu_limit is not None:
            candidates.append(max(1, math.floor(self.cgroup_cpu_limit)))
        return min(candidates) if candidates else None

    def apply_memory(self, memory_bytes: int | None) -> int | None:
        candidates = [value for value in (memory_bytes, self.cgroup_memory_limit_bytes) if value]
        return min(candidates) if candidates else None

    def as_dict(self) -> dict[str, Any]:
        return {
            "affinity_cpus": self.affinity_cpus,
            "cgroup_cpu_limit": self.cgroup_cpu_limit,
            "cgroup_memory_limit_gib": _gib(self.cgroup_memory_limit_bytes),
        }


@dataclass(frozen=True)
class ProcessTuning:
    selector: str
    memory_gib: int
    reason: str

    def as_dict(self) -> dict[str, Any]:
        return {"selector": self.selector, "memory_gib": self.memory_gib, "cpus": "unchanged (nf-core default)", "reason": self.reason}


@dataclass(frozen=True)
class ResourcePolicy:
    """The auditable decision recorded in every run's provenance."""

    cpu_mode: str
    memory_mode: str
    requested_cpus: int
    requested_memory_gib: int
    reserve_cpus: int | None
    reserve_memory_gib: int | None
    detected: DetectedResources
    fallback: bool
    process_tuning: tuple[ProcessTuning, ...] = ()
    notes: tuple[str, ...] = field(default=())

    def with_process_tuning(self, tuning: tuple[ProcessTuning, ...], notes: tuple[str, ...] = ()) -> "ResourcePolicy":
        return ResourcePolicy(
            self.cpu_mode, self.memory_mode, self.requested_cpus, self.requested_memory_gib,
            self.reserve_cpus, self.reserve_memory_gib, self.detected, self.fallback,
            tuning, self.notes + notes,
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": "1.0",
            "cpu_mode": self.cpu_mode,
            "memory_mode": self.memory_mode,
            "detected": self.detected.as_dict(),
            "os_reserve": {"cpus": self.reserve_cpus, "memory_gib": self.reserve_memory_gib},
            "selected": {"cpus": self.requested_cpus, "memory_gib": self.requested_memory_gib},
            "fallback": self.fallback,
            "process_tuning": [item.as_dict() for item in self.process_tuning],
            "notes": list(self.notes + self.detected.warnings),
        }


def _gib(value: int | None) -> float | None:
    return None if value is None else round(value / GIB, 1)


def _read_text(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError):
        return None


def _cgroup_v2_dirs(root: Path) -> list[Path]:
    """This process's cgroup v2 directory and its ancestors, innermost first."""

    base = root / "sys" / "fs" / "cgroup"
    relative = "/"
    for line in (_read_text(root / "proc" / "self" / "cgroup") or "").splitlines():
        if line.startswith("0::"):
            relative = line[3:] or "/"
            break
    current = base / relative.lstrip("/")
    directories: list[Path] = []
    while True:
        directories.append(current)
        if current == base or base not in current.parents:
            break
        current = current.parent
    if base not in directories:
        directories.append(base)
    return directories


def _cgroup_limits(root: Path) -> tuple[float | None, int | None]:
    cpu_limit: float | None = None
    memory_limit: int | None = None
    for directory in _cgroup_v2_dirs(root):
        cpu = _read_text(directory / "cpu.max")
        if cpu:
            quota, _, period = cpu.partition(" ")
            if quota != "max" and quota.isdigit() and period.isdigit() and int(period) > 0:
                value = int(quota) / int(period)
                cpu_limit = value if cpu_limit is None else min(cpu_limit, value)
        memory = _read_text(directory / "memory.max")
        if memory and memory.isdigit() and int(memory) < _UNLIMITED_BYTES:
            memory_limit = int(memory) if memory_limit is None else min(memory_limit, int(memory))
    v1 = root / "sys" / "fs" / "cgroup"
    quota, period = _read_text(v1 / "cpu" / "cpu.cfs_quota_us"), _read_text(v1 / "cpu" / "cpu.cfs_period_us")
    if quota and period and quota.lstrip("-").isdigit() and period.isdigit() and int(quota) > 0 and int(period) > 0:
        value = int(quota) / int(period)
        cpu_limit = value if cpu_limit is None else min(cpu_limit, value)
    memory = _read_text(v1 / "memory" / "memory.limit_in_bytes")
    if memory and memory.isdigit() and int(memory) < _UNLIMITED_BYTES:
        memory_limit = int(memory) if memory_limit is None else min(memory_limit, int(memory))
    return cpu_limit, memory_limit


def detect_process_limits(root: Path = Path("/")) -> ProcessLimits:
    """CPU affinity and cgroup v1/v2 limits for this process; never raises."""

    affinity: int | None = None
    try:
        affinity = len(os.sched_getaffinity(0)) or None  # type: ignore[attr-defined]
    except (AttributeError, OSError):
        affinity = None
    try:
        cgroup_cpu, cgroup_memory = _cgroup_limits(root)
    except Exception:  # detection must never block a run
        cgroup_cpu, cgroup_memory = None, None
    return ProcessLimits(affinity, cgroup_cpu, cgroup_memory)


def cpu_reserve(cpus: int) -> int:
    return max(CPU_RESERVE_MINIMUM, math.ceil(cpus * CPU_RESERVE_FRACTION))


def memory_reserve_gib(memory_gib: float) -> int:
    return max(MEMORY_RESERVE_MINIMUM_GIB, math.ceil(memory_gib * MEMORY_RESERVE_FRACTION))


def resolve_policy(configured_cpus: int | str, configured_memory_gib: int | str, detected: DetectedResources) -> ResourcePolicy:
    """Select the aggregate local ceiling.

    Explicit integers from ``project.yaml`` are returned unchanged; the
    existing effective-budget step may still clamp them to physical capacity
    with a visible warning.  ``auto`` uses detected capacity minus an OS
    reserve, never below the 8 CPU / 12 GiB contract floor unless the machine
    itself is smaller.  Undetectable capacity falls back to 8 CPUs / 12 GiB.
    """

    fallback = False
    notes: list[str] = []
    reserve_cpus: int | None = None
    reserve_memory: int | None = None
    if configured_cpus == AUTO:
        usable = detected.usable_cpus
        if usable is None:
            cpus, fallback = FALLBACK_CPUS, True
            notes.append(f"CPU capacity undetectable; using the fallback of {FALLBACK_CPUS} CPUs.")
        else:
            reserve_cpus = cpu_reserve(usable)
            cpus = min(usable, max(FLOOR_CPUS, usable - reserve_cpus))
            if cpus > usable - reserve_cpus:
                notes.append(f"CPU reserve reduced to {usable - cpus} to satisfy the {FLOOR_CPUS}-CPU local contract floor.")
    else:
        cpus = int(configured_cpus)
    if configured_memory_gib == AUTO:
        usable_bytes = detected.usable_memory_bytes
        if usable_bytes is None:
            memory, fallback = FALLBACK_MEMORY_GIB, True
            notes.append(f"Memory capacity undetectable; using the fallback of {FALLBACK_MEMORY_GIB} GiB.")
        else:
            usable_gib = usable_bytes / GIB
            reserve_memory = memory_reserve_gib(usable_gib)
            memory = min(math.floor(usable_gib), max(FLOOR_MEMORY_GIB, math.floor(usable_gib - reserve_memory)))
            memory = max(1, memory)
            if memory > usable_gib - reserve_memory:
                notes.append(f"Memory reserve reduced to satisfy the {FLOOR_MEMORY_GIB}-GiB local contract floor.")
            available = detected.available_memory_bytes
            if available is not None and available / GIB < memory:
                notes.append(
                    f"Only {available / GIB:.1f} GiB is currently available (selected {memory} GiB); "
                    "other workloads may compete for memory."
                )
    else:
        memory = int(configured_memory_gib)
    return ResourcePolicy(
        cpu_mode=AUTO if configured_cpus == AUTO else "explicit",
        memory_mode=AUTO if configured_memory_gib == AUTO else "explicit",
        requested_cpus=cpus,
        requested_memory_gib=memory,
        reserve_cpus=reserve_cpus,
        reserve_memory_gib=reserve_memory,
        detected=detected,
        fallback=fallback,
        notes=tuple(notes),
    )


def describe_policy(policy: ResourcePolicy) -> str:
    """One line for doctor/status: how each ceiling was chosen."""

    def part(label: str, mode: str, value: int, usable: str, reserve: int | None, unit: str) -> str:
        if mode != AUTO:
            return f"{label} {value}{unit} explicit in project.yaml"
        if reserve is None:
            return f"{label} {value}{unit} auto fallback (capacity undetectable)"
        return f"{label} {value}{unit} auto (usable {usable}{unit}, OS reserve {reserve}{unit})"

    usable_memory = policy.detected.usable_memory_bytes
    return "; ".join((
        part("CPUs", policy.cpu_mode, policy.requested_cpus, str(policy.detected.usable_cpus or "?"), policy.reserve_cpus, ""),
        part("memory", policy.memory_mode, policy.requested_memory_gib,
             f"{usable_memory / GIB:.0f}" if usable_memory else "?", policy.reserve_memory_gib, " GiB"),
    ))


def directory_size_bytes(path: Path) -> int:
    total = 0
    for item in path.rglob("*"):
        if item.is_file() and not item.is_symlink():
            total += item.stat().st_size
    return total


def salmon_quant_tuning(index: Path | None) -> tuple[tuple[ProcessTuning, ...], tuple[str, ...]]:
    """Size SALMON_QUANT's memory request from a prebuilt index, or leave nf-core's default."""

    if index is None:
        return (), ("Salmon quantification keeps nf-core's default request: no prebuilt Salmon index to size it from.",)
    try:
        size = directory_size_bytes(index) if index.is_dir() else 0
    except OSError as exc:
        return (), (f"Salmon quantification keeps nf-core's default request: index size unavailable ({exc}).",)
    if size <= 0:
        return (), ("Salmon quantification keeps nf-core's default request: index size unavailable.",)
    derived = math.ceil(size / GIB * SALMON_INDEX_MEMORY_FACTOR + SALMON_INDEX_MEMORY_OVERHEAD_GIB)
    memory = min(SALMON_QUANT_NFCORE_MEMORY_GIB, derived)
    return (ProcessTuning(
        SALMON_QUANT_SELECTOR, memory,
        f"index {size / GIB:.1f} GiB x {SALMON_INDEX_MEMORY_FACTOR} + {SALMON_INDEX_MEMORY_OVERHEAD_GIB} GiB, "
        f"capped at nf-core's {SALMON_QUANT_NFCORE_MEMORY_GIB} GiB",
    ),), ()


def render_process_tuning_config(tuning: tuple[ProcessTuning, ...]) -> str:
    """Scheduling-only per-process requests for the nf-core upstream workflow."""

    lines = [
        "// nf-rna scheduling requests only: memory requests let independent tasks run concurrently.",
        "// CPUs (tool thread counts) and every tool argument are unchanged.",
        "process {",
    ]
    for item in tuning:
        lines.append(f"  withName: {item.selector!r} {{ memory = '{item.memory_gib}.GB' }}")
    lines.append("}")
    return "\n".join(lines) + "\n"
