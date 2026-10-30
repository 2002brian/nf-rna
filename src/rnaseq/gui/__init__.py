"""Optional browser GUI (the ``gui`` extra; started with ``rnaseq gui``).

Only :mod:`rnaseq.gui.app` and :mod:`rnaseq.gui.pages` import Streamlit; the
service modules here are plain Python so they stay importable and testable
without the extra.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

GUI_ADDRESS = "127.0.0.1"
DEFAULT_PORT = 8501
INSTALL_HINT = "pip install -e '.[gui]'  (from a source checkout)  or  python -m pip install \"nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@vX.Y.Z\""


def streamlit_available() -> bool:
    return importlib.util.find_spec("streamlit") is not None


def streamlit_command(port: int) -> list[str]:
    """The Streamlit server command; it listens on the loopback interface only."""

    return [
        sys.executable, "-m", "streamlit", "run", str(Path(__file__).with_name("app.py")),
        "--server.address", GUI_ADDRESS, "--server.port", str(port), "--server.headless", "true",
        "--browser.gatherUsageStats", "false", "--server.fileWatcherType", "none", "--client.toolbarMode", "viewer",
    ]


def launch(port: int) -> None:
    """Replace this process with the Streamlit server (so a service manager supervises it directly)."""

    command = streamlit_command(port)
    os.execv(command[0], command)
