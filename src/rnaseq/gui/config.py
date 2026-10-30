"""GUI settings (``~/.config/nf-rna/gui.yaml``) and the allowed-roots path guard."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import yaml

# Written to gui.yaml on first use; ``~`` is the user's home directory.
DEFAULT_ALLOWED_ROOTS = ("/media/duobook02/brian/nf-rna", "~/nf-rna-work/testdata")
DEFAULT_CONDARC = "~/nf-rna-work/conda/condarc.yml"
# Roots below these prefixes live on removable or network storage, so an
# unmounted disk must not be mistaken for an empty local directory.
MOUNT_PREFIXES = ("/media/", "/mnt/", "/run/media/")


class PathNotAllowed(ValueError):
    """A path outside the configured allowed roots (after resolving symlinks)."""


def config_path() -> Path:
    if os.environ.get("NF_RNA_GUI_CONFIG"):
        return Path(os.environ["NF_RNA_GUI_CONFIG"]).expanduser()
    base = Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config")
    return base / "nf-rna" / "gui.yaml"


def gui_state_dir() -> Path:
    """GUI launch records, the active-run lock and launcher logs (never analysis state)."""

    base = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
    return base / "nf-rna" / "gui"


@dataclass(frozen=True)
class RootStatus:
    path: Path
    available: bool
    message: str


@dataclass(frozen=True)
class GuiConfig:
    allowed_roots: tuple[Path, ...]
    project_parent: Path
    condarc: Path | None
    log_tail_lines: int = 200
    launcher: str = "auto"
    # Larger files and delivery packages are offered as a server path and scp command.
    download_limit_mb: int = 500
    source: Path | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)

    def root_statuses(self) -> list[RootStatus]:
        return [root_status(root) for root in self.allowed_roots]

    def available_roots(self) -> list[Path]:
        return [status.path for status in self.root_statuses() if status.available]

    def resolve_allowed(self, value: str | Path) -> Path:
        """Resolve ``value`` and require it to lie inside an available allowed root."""

        text = str(value).strip()
        if not text:
            raise PathNotAllowed("No path given.")
        path = Path(text).expanduser().resolve()
        for root in self.available_roots():
            resolved_root = root.resolve()
            if path == resolved_root or resolved_root in path.parents:
                return path
        raise PathNotAllowed(f"{path} is outside the allowed roots: {', '.join(str(root) for root in self.allowed_roots)}.")


def _mount_point(path: Path) -> Path:
    current = path
    while not os.path.ismount(current) and current != current.parent:
        current = current.parent
    return current


def root_status(root: Path) -> RootStatus:
    if any(str(root).startswith(prefix) for prefix in MOUNT_PREFIXES):
        mount = _mount_point(root if root.exists() else root.parent if root.parent.exists() else Path("/"))
        if mount == Path("/"):
            return RootStatus(root, False, f"not mounted: {root} would be on the root filesystem; mount the data disk first.")
        if not root.is_dir():
            return RootStatus(root, False, f"not found on the mounted filesystem {mount}.")
        return RootStatus(root, True, f"mounted ({mount})")
    if not root.is_dir():
        return RootStatus(root, False, "does not exist.")
    return RootStatus(root, True, "available")


def load_config(path: Path | None = None) -> GuiConfig:
    """Read gui.yaml; a missing file is written once with the documented defaults."""

    path = path or config_path()
    raw: dict = {}
    notes: list[str] = []
    if path.is_file():
        loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
        raw = loaded if isinstance(loaded, dict) else {}
    else:
        condarc = str(Path(DEFAULT_CONDARC).expanduser()) if Path(DEFAULT_CONDARC).expanduser().is_file() else None
        raw = {"allowed_roots": [str(Path(root).expanduser()) for root in DEFAULT_ALLOWED_ROOTS], "project_parent": f"{DEFAULT_ALLOWED_ROOTS[0]}/projects", "condarc": condarc, "log_tail_lines": 200, "launcher": "auto", "download_limit_mb": 500}
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
            notes.append(f"Wrote default GUI settings to {path}.")
        except OSError as exc:
            notes.append(f"Using default GUI settings; could not write {path}: {exc}")
    roots = tuple(Path(str(item)).expanduser() for item in raw.get("allowed_roots") or DEFAULT_ALLOWED_ROOTS)
    parent = Path(str(raw.get("project_parent") or roots[0] / "projects")).expanduser()
    condarc = Path(str(raw["condarc"])).expanduser() if raw.get("condarc") else None
    launcher = str(raw.get("launcher") or "auto")
    if launcher not in {"auto", "systemd", "process"}:
        raise ValueError(f"{path}: launcher must be auto, systemd or process (got {launcher!r}).")
    return GuiConfig(roots, parent, condarc, int(raw.get("log_tail_lines") or 200), launcher, int(raw.get("download_limit_mb") or 500), path, tuple(notes))


def ensure_data_layout(config: GuiConfig) -> list[str]:
    """Create ``projects/`` and ``fastq/`` beside each other in the data root, only when it is mounted."""

    data_root = config.project_parent.parent
    status = next((item for item in config.root_statuses() if item.path == data_root), root_status(data_root))
    if not status.available:
        return [f"{data_root}: {status.message} Project and FASTQ folders were not created."]
    created = []
    for directory in (config.project_parent, data_root / "fastq"):
        if not directory.is_dir():
            directory.mkdir(parents=True, exist_ok=True)
            created.append(f"Created {directory}.")
    return created
