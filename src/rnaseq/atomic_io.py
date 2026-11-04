"""Crash-safe replacement of small control-plane files.

Run state, frozen contracts and manifests are polled by ``rnaseq status``, the
GUI and later runs.  A plain ``write_text`` truncates the file first, so an
interruption (SIGKILL, OOM, power loss) or a concurrent reader can observe an
empty or partial file.  Writing a sibling temporary file, flushing it to disk
and renaming it over the target makes every reader see either the previous or
the new complete content.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path


def _umask_mode() -> int:
    mask = os.umask(0)
    os.umask(mask)
    return 0o666 & ~mask


# Computed once: os.umask() can only be read by setting it, which is not
# thread-safe to repeat while other threads create files.
_DEFAULT_MODE = _umask_mode()


def atomic_write_text(path: Path, content: str) -> None:
    """Atomically replace ``path`` with UTF-8 ``content`` using LF newlines.

    The file keeps the permission bits a plain ``write_text`` would have given
    it (the existing file's mode, else the process umask), not ``mkstemp``'s
    private 0600.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        mode = path.stat().st_mode & 0o7777
    except OSError:
        mode = _DEFAULT_MODE
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    try:  # Persist the rename itself; not every platform can fsync a directory.
        directory = os.open(path.parent, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(directory)
    except OSError:
        pass
    finally:
        os.close(directory)
