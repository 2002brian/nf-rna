"""Provision the immutable native runtime for first-party downstream tasks.

This module intentionally owns only the nf-rna downstream Conda environment.
Upstream nf-core and featureCounts containers remain outside this contract.

The nf-rna wheel installed into the locked prefix comes from one of two sources:

* a normal (non-editable) installation: the installed distribution is verified
  against its pip ``RECORD`` hashes and re-packed, with the standard library
  only, into a deterministic canonical wheel; no build tool, network, or Git
  checkout is involved;
* a development source checkout (including editable installs): the wheel is
  built from the checkout with ``python -m build`` (the ``dev`` extra).
"""

from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from importlib import metadata, resources
from pathlib import Path
from typing import Any

import yaml

from rnaseq.errors import ExecutionPreflightError


RUNTIME_SCHEMA = "nf-rna.downstream-conda-runtime.v1"
DISTRIBUTION_NAME = "nf-rna"
WHEEL_ORIGIN_SOURCE = "source-checkout-build"
WHEEL_ORIGIN_INSTALLED = "installed-distribution"
# pip writes these at install time; they are not part of the distributed wheel.
INSTALLER_GENERATED = frozenset({"INSTALLER", "REQUESTED", "direct_url.json", "RECORD"})
CANONICAL_WHEEL_TIMESTAMP = (2020, 2, 2, 0, 0, 0)
LOCKS_PACKAGE = "rnaseq.workflows"
REQUIRED_R_PACKAGES = (
    "DESeq2", "tximport", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db",
    "org.Mm.eg.db", "ggplot2", "pheatmap", "yaml", "jsonlite",
)


@dataclass(frozen=True)
class DownstreamRuntime:
    """The verified prefix and immutable identities consumed by Nextflow."""

    prefix: Path
    platform: str
    lock_filename: str
    lock_sha256: str
    wheel_filename: str
    wheel_sha256: str
    nf_rna_version: str
    source_revision: str | None
    r_scripts: tuple[dict[str, str], ...]
    r_scripts_sha256: str
    wheel_origin: str = WHEEL_ORIGIN_SOURCE

    def identity(self) -> dict[str, Any]:
        return {
            "schema": RUNTIME_SCHEMA,
            "kind": "conda",
            "platform": self.platform,
            "prefix": str(self.prefix),
            "lock": {"filename": self.lock_filename, "sha256": self.lock_sha256},
            "wheel": {"filename": self.wheel_filename, "sha256": self.wheel_sha256, "origin": self.wheel_origin},
            "nf_rna_version": self.nf_rna_version,
            "source_revision": self.source_revision,
            "r_scripts": {"inventory": list(self.r_scripts), "sha256": self.r_scripts_sha256},
        }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def native_platform() -> str:
    """Return the explicit lock platform for the native machine."""

    machine = platform.machine().lower()
    system = platform.system().lower()
    if system == "darwin" and machine in {"arm64", "aarch64"}:
        return "osx-arm64"
    if system == "linux" and machine in {"x86_64", "amd64"}:
        return "linux-64"
    raise ExecutionPreflightError(
        f"No native nf-rna downstream Conda lock is available for {system}-{machine}."
    )


def _workflow_path(relative: str) -> Path:
    item = resources.files(LOCKS_PACKAGE).joinpath(relative)
    try:
        path = Path(item)
    except TypeError as exc:  # pragma: no cover - wheels install filesystem resources.
        raise ExecutionPreflightError("Bundled runtime contract is not filesystem-backed.") from exc
    if not path.is_file():
        # Hatch maps workflow/ into wheels.  A source checkout intentionally
        # keeps those assets at repository level for the control plane only;
        # downstream task code and R scripts still come from the installed
        # wheel provisioned below.
        development = Path(__file__).resolve().parents[2] / "workflow" / relative
        if development.is_file():
            return development
        raise ExecutionPreflightError(f"Bundled runtime contract is missing: {relative}.")
    return path


def lock_path_for_platform(target: str | None = None) -> Path:
    target = target or native_platform()
    if target not in {"linux-64", "osx-arm64"}:
        raise ExecutionPreflightError(f"Unsupported downstream Conda platform: {target}.")
    return _workflow_path(f"envs/locks/nf-rna-downstream-{target}.lock.yml")


def _lock_checksum(lock: Path) -> str:
    sums = _workflow_path("envs/locks/SHA256SUMS")
    expected: str | None = None
    for line in sums.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[1] == lock.name:
            expected = fields[0]
            break
    if expected is None:
        raise ExecutionPreflightError(f"No recorded SHA-256 for downstream runtime lock {lock.name}.")
    observed = _sha256(lock)
    if observed != expected:
        raise ExecutionPreflightError(
            f"Downstream runtime lock checksum mismatch for {lock.name}: expected {expected}, observed {observed}."
        )
    return observed


def _runtime_root() -> Path:
    configured = os.environ.get("RNASEQ_RUNTIME_ROOT")
    if configured:
        return Path(configured).expanduser().resolve()
    if platform.system() == "Darwin":
        return Path.home() / "Library" / "Caches" / "nf-rna" / "runtime"
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "nf-rna" / "runtime"


def _source_checkout() -> Path | None:
    """Return the repository root only when this module runs from its ``src/`` tree."""

    module = Path(__file__).resolve()
    root = module.parents[2]
    if module.parents[1].name == "src" and (root / "pyproject.toml").is_file():
        return root
    return None


def _source_root() -> Path:
    root = _source_checkout()
    if root is None:
        raise ExecutionPreflightError("This nf-rna invocation is not a source checkout.")
    return root


def _installed_distribution() -> metadata.Distribution:
    """Return the installed distribution that owns the running ``rnaseq`` package."""

    import rnaseq

    package = Path(rnaseq.__file__).resolve()
    for distribution in metadata.distributions(name=DISTRIBUTION_NAME):
        if Path(distribution.locate_file("rnaseq/__init__.py")).resolve() == package:
            return distribution
    raise ExecutionPreflightError(
        f"The running nf-rna package ({package.parent}) is neither a source checkout nor owned by an "
        "installed nf-rna distribution; install nf-rna with pip (for example 'pip install nf-rna')."
    )


def _record_sha256(data: bytes) -> str:
    return "sha256=" + base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode("ascii")


def _installed_wheel_members(distribution: metadata.Distribution) -> tuple[str, dict[str, bytes]]:
    """Return the verified distributed files of an installed nf-rna.

    Every packaged file must still match the SHA-256 pip recorded at install
    time, so a locally modified installation is refused rather than provisioned.
    """

    record = distribution.read_text("RECORD")
    if not record:
        raise ExecutionPreflightError("The installed nf-rna distribution has no RECORD; reinstall nf-rna with pip.")
    rows = list(csv.reader(io.StringIO(record)))
    info_dirs = {row[0].split("/", 1)[0] for row in rows if row and row[0].endswith(".dist-info/METADATA") and row[0].count("/") == 1}
    if len(info_dirs) != 1:
        raise ExecutionPreflightError("The installed nf-rna RECORD does not identify exactly one .dist-info directory.")
    info_dir = info_dirs.pop()
    members: dict[str, bytes] = {}
    for row in rows:
        if len(row) != 3:
            raise ExecutionPreflightError(f"The installed nf-rna RECORD has a malformed row: {row!r}.")
        path, digest, _size = row
        top, _, name = path.partition("/")
        if top not in {"rnaseq", info_dir} or "__pycache__/" in path or path.endswith(".pyc"):
            continue
        if top == info_dir and name in INSTALLER_GENERATED:
            continue
        if not digest.startswith("sha256="):
            raise ExecutionPreflightError(f"The installed nf-rna RECORD has no SHA-256 for {path}; reinstall nf-rna.")
        try:
            data = Path(distribution.locate_file(path)).read_bytes()
        except OSError as exc:
            raise ExecutionPreflightError(f"Installed nf-rna file is missing: {path}; reinstall nf-rna.") from exc
        if _record_sha256(data) != digest:
            raise ExecutionPreflightError(f"Installed nf-rna file was modified after installation: {path}; reinstall nf-rna.")
        members[path] = data
    required = ("rnaseq/workflow_support.py", "rnaseq/workflows/main.nf", "rnaseq/workflows/envs/locks/SHA256SUMS", f"{info_dir}/METADATA", f"{info_dir}/WHEEL")
    missing = [path for path in required if path not in members]
    if missing or not any(path.startswith("rnaseq/r/") and path.endswith(".R") for path in members):
        raise ExecutionPreflightError("The installed nf-rna distribution is incomplete: " + ", ".join(missing or ["rnaseq/r/*.R"]))
    return info_dir, members


def _canonical_wheel_bytes(info_dir: str, members: dict[str, bytes]) -> bytes:
    """Deterministic, uncompressed wheel: identical content always yields identical bytes."""

    ordered = sorted(path for path in members if not path.startswith(info_dir + "/"))
    ordered += sorted(path for path in members if path.startswith(info_dir + "/"))
    record = "".join(f"{path},{_record_sha256(members[path])},{len(members[path])}\n" for path in ordered)
    entries = [(path, members[path]) for path in ordered] + [(f"{info_dir}/RECORD", (record + f"{info_dir}/RECORD,,\n").encode())]
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for path, data in entries:
            info = zipfile.ZipInfo(path, CANONICAL_WHEEL_TIMESTAMP)
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)
    return buffer.getvalue()


def _installed_wheel_filename(info_dir: str, members: dict[str, bytes]) -> str:
    wheel_metadata = members[f"{info_dir}/WHEEL"].decode("utf-8")
    tags = [line.split(":", 1)[1].strip() for line in wheel_metadata.splitlines() if line.startswith("Tag:")]
    if len(tags) != 1:
        raise ExecutionPreflightError(f"The installed nf-rna WHEEL metadata must declare exactly one tag; found {tags}.")
    return f"{info_dir.removesuffix('.dist-info')}-{tags[0]}.whl"


def _installed_revision(distribution: metadata.Distribution) -> str | None:
    """Use the VCS commit pip recorded for a Git-URL install; releases carry none."""

    try:
        direct = json.loads(distribution.read_text("direct_url.json") or "{}")
    except json.JSONDecodeError:
        return None
    commit = direct.get("vcs_info", {}).get("commit_id") if isinstance(direct.get("vcs_info"), dict) else None
    return commit if isinstance(commit, str) and commit else None


def _materialize_installed_wheel(cache: Path) -> tuple[Path, str | None]:
    """Re-pack the verified installed distribution into a content-addressed canonical wheel."""

    distribution = _installed_distribution()
    info_dir, members = _installed_wheel_members(distribution)
    data = _canonical_wheel_bytes(info_dir, members)
    digest = hashlib.sha256(data).hexdigest()
    target = cache / "wheels" / "installed" / digest / _installed_wheel_filename(info_dir, members)
    if not target.is_file() or _sha256(target) != digest:
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as handle:
            handle.write(data)
        os.replace(handle.name, target)
    return target.resolve(), _installed_revision(distribution)


def _runtime_wheel(cache: Path) -> tuple[Path, str | None, str]:
    """Select the nf-rna wheel for the downstream prefix without needing a checkout in production."""

    if _source_checkout() is not None:
        wheel, revision = _build_wheel(cache)
        return wheel, revision, WHEEL_ORIGIN_SOURCE
    wheel, revision = _materialize_installed_wheel(cache)
    return wheel, revision, WHEEL_ORIGIN_INSTALLED


def wheel_source_status() -> tuple[bool, str]:
    """Non-mutating doctor view of where the downstream nf-rna wheel will come from."""

    source = _source_checkout()
    if source is not None:
        import importlib.util

        if importlib.util.find_spec("build") is None:
            return False, (
                f"source checkout {source} builds the downstream wheel with '{sys.executable} -m build', which is "
                "not importable; install the 'dev' extra, or install nf-rna normally (non-editable) for production."
            )
        return True, f"source checkout {source}; development wheel built with '{sys.executable} -m build'."
    try:
        distribution = _installed_distribution()
        info_dir, members = _installed_wheel_members(distribution)
    except ExecutionPreflightError as exc:
        return False, str(exc)
    return True, (
        f"installed distribution {info_dir.removesuffix('.dist-info')}; {len(members)} files verified against RECORD; "
        "canonical wheel is re-packed without build tools."
    )


def _source_revision(source: Path) -> str | None:
    if not (source / ".git").exists():
        return None
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=source, capture_output=True, text=True, check=False
    )
    revision = result.stdout.strip() if result.returncode == 0 and result.stdout.strip() else None
    if revision is None:
        return None
    status = subprocess.run(
        ["git", "status", "--porcelain"], cwd=source, capture_output=True, text=True, check=False
    )
    return revision + "+dirty" if status.returncode == 0 and status.stdout.strip() else revision


def _build_wheel(cache: Path) -> tuple[Path, str | None]:
    source = _source_root()
    revision = _source_revision(source)
    wheel_dir = cache / "wheels"
    wheel_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(wheel_dir), str(source)],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        raise ExecutionPreflightError(
            f"Could not build the nf-rna downstream wheel from source checkout {source}: {detail} "
            "(development mode needs the 'dev' extra; production installs do not build wheels)."
        )
    wheels = sorted(wheel_dir.glob("nf_rna-*.whl"), key=lambda item: item.stat().st_mtime_ns)
    if not wheels:
        raise ExecutionPreflightError("Wheel build completed without an nf-rna wheel artifact.")
    return wheels[-1].resolve(), revision


def _r_script_inventory(prefix: Path) -> tuple[tuple[dict[str, str], ...], str]:
    python = prefix / "bin" / "python"
    located = subprocess.run(
        [str(python), "-c", "from importlib.resources import files; print(files('rnaseq').joinpath('r'))"],
        capture_output=True, text=True, check=False,
    )
    candidates = [Path(line.strip()).resolve() for line in located.stdout.splitlines() if line.strip()]
    if located.returncode != 0 or len(candidates) != 1 or prefix.resolve() not in candidates[0].parents:
        raise ExecutionPreflightError("Installed nf-rna wheel has no unique bundled R-script directory.")
    scripts = tuple(
        {"name": path.name, "sha256": _sha256(path)}
        for path in sorted(candidates[0].glob("*.R"))
    )
    if not scripts:
        raise ExecutionPreflightError("Installed nf-rna wheel does not contain bundled R scripts.")
    digest = hashlib.sha256(json.dumps(scripts, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return scripts, digest


def _probe(prefix: Path, expected_wheel: Path) -> tuple[str, tuple[dict[str, str], ...], str]:
    python = prefix / "bin" / "python"
    rscript = prefix / "bin" / "Rscript"
    if not python.is_file() or not rscript.is_file():
        raise ExecutionPreflightError("Locked downstream Conda prefix is missing Python or Rscript.")
    script = """
import json
import pathlib
import rnaseq
import rnaseq.workflow_support
root = pathlib.Path(rnaseq.__file__).resolve()
prefix = pathlib.Path(__import__('sys').prefix).resolve()
if prefix not in root.parents:
    raise SystemExit(f'nf-rna resolved outside Conda prefix: {root}')
print(json.dumps({'version': rnaseq.__version__, 'file': str(root)}))
"""
    result = subprocess.run([str(python), "-c", script], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise ExecutionPreflightError(f"Installed downstream wheel failed Python probe: {result.stderr.strip()}")
    try:
        identity = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise ExecutionPreflightError("Installed downstream wheel Python probe returned invalid data.") from exc
    r_code = "missing <- c(" + ",".join(json.dumps(name) for name in REQUIRED_R_PACKAGES) + "); bad <- missing[!vapply(missing, requireNamespace, logical(1), quietly=TRUE)]; if (length(bad)) stop(paste(bad, collapse=','))"
    r_result = subprocess.run([str(rscript), "-e", r_code], capture_output=True, text=True, check=False)
    if r_result.returncode != 0:
        raise ExecutionPreflightError(f"Locked downstream R package probe failed: {r_result.stderr.strip()}")
    scripts, digest = _r_script_inventory(prefix)
    installed_wheel = prefix / "runtime" / "nf-rna-wheel.json"
    if installed_wheel.is_file():
        recorded = json.loads(installed_wheel.read_text(encoding="utf-8"))
        expected_sha = _sha256(expected_wheel)
        if recorded.get("sha256") != expected_sha:
            raise ExecutionPreflightError(
                "Cached downstream runtime wheel identity does not match the requested wheel: "
                f"prefix {prefix} holds {recorded.get('filename')} sha256={recorded.get('sha256')}, "
                f"requested {expected_wheel.name} sha256={expected_sha}. The prefix is never updated in place; "
                "remove it manually or set RNASEQ_RUNTIME_ROOT to a separate runtime root."
            )
    return str(identity["version"]), scripts, digest


def downstream_runtime_preflight() -> None:
    """Verify immutable inputs without creating a Conda environment."""

    lock = lock_path_for_platform()
    _lock_checksum(lock)
    if not shutil_which("conda"):
        raise ExecutionPreflightError("Conda is required to provision the locked downstream runtime.")


def shutil_which(name: str) -> str | None:
    # Kept tiny to make preflight deterministic and easy to replace in tests.
    import shutil
    return shutil.which(name)


def ensure_downstream_runtime() -> DownstreamRuntime:
    """Create or verify the cached non-editable Conda runtime after confirmation."""

    downstream_runtime_preflight()
    target = native_platform()
    lock = lock_path_for_platform(target)
    lock_sha = _lock_checksum(lock)
    root = _runtime_root().resolve()
    prefix = root / "prefixes" / f"nf-rna-downstream-{target}-{lock_sha[:16]}"
    wheel, revision, origin = _runtime_wheel(root)
    wheel_sha = _sha256(wheel)
    marker = prefix / "runtime" / "nf-rna-wheel.json"
    if not prefix.is_dir():
        env = {**os.environ, "CONDA_NO_PLUGINS": "true", "CONDA_SOLVER": "classic"}
        result = subprocess.run(
            ["conda", "env", "create", "--yes", "--prefix", str(prefix), "--file", str(lock)],
            capture_output=True, text=True, check=False, env=env,
        )
        if result.returncode != 0:
            detail = (result.stderr or result.stdout).strip()
            raise ExecutionPreflightError(f"Could not create locked downstream Conda environment: {detail}")
        install = subprocess.run(
            [str(prefix / "bin" / "python"), "-m", "pip", "install", "--no-deps", str(wheel)],
            capture_output=True, text=True, check=False,
        )
        if install.returncode != 0:
            raise ExecutionPreflightError(f"Could not install non-editable nf-rna wheel: {(install.stderr or install.stdout).strip()}")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(json.dumps({"filename": wheel.name, "sha256": wheel_sha, "origin": origin}, indent=2) + "\n", encoding="utf-8")
    elif not marker.is_file():
        raise ExecutionPreflightError(
            f"Existing downstream runtime prefix lacks its wheel identity marker: {prefix}. Remove it manually and retry."
        )
    version, scripts, scripts_sha = _probe(prefix, wheel)
    return DownstreamRuntime(prefix, target, lock.name, lock_sha, wheel.name, wheel_sha, version, revision, scripts, scripts_sha, origin)
