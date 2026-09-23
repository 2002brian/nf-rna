"""Production packaging contract: a normal installed nf-rna provisions its downstream runtime.

The synthetic installation mirrors what ``pip install nf_rna-*.whl`` writes:
the package tree (with workflow/ mapped to rnaseq/workflows/), a .dist-info
directory, and a RECORD that also lists pip-generated files, bytecode, and the
console script outside site-packages.
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile
from importlib import metadata
from pathlib import Path

import pytest

from rnaseq import downstream_runtime
from rnaseq.errors import ExecutionPreflightError


ROOT = Path(__file__).parents[1]
INFO = "nf_rna-1.2.0.dist-info"


def _record_hash(data: bytes) -> str:
    return "sha256=" + base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()


def _write_record(site: Path) -> None:
    rows = []
    for path in sorted(site.rglob("*")):
        relative = path.relative_to(site).as_posix()
        if path.is_dir() or relative == f"{INFO}/RECORD":
            continue
        rows.append(f"{relative},," if "__pycache__" in relative else f"{relative},{_record_hash(path.read_bytes())},{path.stat().st_size}")
    rows += ["../../../bin/rnaseq,sha256=AAAA,10", f"{INFO}/RECORD,,"]
    (site / INFO / "RECORD").write_text("\n".join(rows) + "\n", encoding="utf-8")


def _install(site: Path, *, direct_url: dict | None = None) -> metadata.Distribution:
    """Lay out a non-editable nf-rna installation below ``site``."""

    ignore = shutil.ignore_patterns("__pycache__")
    shutil.copytree(ROOT / "src" / "rnaseq", site / "rnaseq", ignore=ignore)
    shutil.copytree(ROOT / "workflow", site / "rnaseq" / "workflows", ignore=ignore, dirs_exist_ok=True)
    (site / "rnaseq" / "__pycache__").mkdir()
    (site / "rnaseq" / "__pycache__" / "cli.cpython-311.pyc").write_bytes(b"bytecode")
    info = site / INFO
    (info / "licenses").mkdir(parents=True)
    (info / "METADATA").write_text("Metadata-Version: 2.4\nName: nf-rna\nVersion: 1.2.0\n", encoding="utf-8")
    (info / "WHEEL").write_text("Wheel-Version: 1.0\nGenerator: hatchling\nRoot-Is-Purelib: true\nTag: py3-none-any\n", encoding="utf-8")
    (info / "entry_points.txt").write_text("[console_scripts]\nrnaseq = rnaseq.cli:app\n", encoding="utf-8")
    (info / "licenses" / "LICENSE").write_bytes((ROOT / "LICENSE").read_bytes())
    (info / "INSTALLER").write_text("pip\n", encoding="utf-8")
    (info / "REQUESTED").write_text("", encoding="utf-8")
    (info / "direct_url.json").write_text(json.dumps(direct_url or {"archive_info": {}, "url": "file:///release.whl"}), encoding="utf-8")
    _write_record(site)
    return metadata.PathDistribution(info)


@pytest.fixture
def installed(monkeypatch, tmp_path):
    """Run the provisioning code as if nf-rna were a normal non-editable install."""

    distribution = _install(tmp_path / "site")
    monkeypatch.setattr(downstream_runtime, "_source_checkout", lambda: None)
    monkeypatch.setattr(downstream_runtime, "_installed_distribution", lambda: distribution)
    monkeypatch.setattr(downstream_runtime, "_build_wheel", lambda *_args: pytest.fail("installed mode built a wheel"))
    return tmp_path / "site"


class FakeConda:
    """Creates a prefix and installs the wheel by extraction; records every command."""

    def __init__(self, monkeypatch):
        self.calls: list[list[str]] = []
        monkeypatch.setattr(downstream_runtime, "shutil_which", lambda name: f"/opt/conda/bin/{name}")
        monkeypatch.setattr(downstream_runtime.subprocess, "run", self)

    def __call__(self, arguments, **_kwargs):
        arguments = [str(item) for item in arguments]
        self.calls.append(arguments)
        assert not any("docker" in item for item in arguments), arguments
        assert "build" not in arguments, "installed mode must not invoke python -m build"
        done = subprocess.CompletedProcess(arguments, 0, "", "")
        if arguments[:3] == ["conda", "env", "create"]:
            prefix = Path(arguments[arguments.index("--prefix") + 1])
            (prefix / "bin").mkdir(parents=True)
            (prefix / "bin" / "python").write_text("", encoding="utf-8")
            (prefix / "bin" / "Rscript").write_text("", encoding="utf-8")
            return done
        prefix = Path(arguments[0]).parents[1]
        site = prefix / "lib" / "site-packages"
        if arguments[1:5] == ["-m", "pip", "install", "--no-deps"]:
            with zipfile.ZipFile(arguments[5]) as wheel:
                wheel.extractall(site)
            return done
        if arguments[1] == "-c" and "print(files('rnaseq').joinpath('r'))" in arguments[2]:
            return subprocess.CompletedProcess(arguments, 0, f"{site / 'rnaseq' / 'r'}\n", "")
        if arguments[1] == "-c":
            return subprocess.CompletedProcess(arguments, 0, json.dumps({"version": "1.2.0", "file": str(site / "rnaseq" / "__init__.py")}), "")
        if arguments[1] == "-e":
            return done
        raise AssertionError(f"unexpected provisioning command: {arguments}")

    def created(self) -> int:
        return sum(call[:3] == ["conda", "env", "create"] for call in self.calls)


def test_installed_distribution_repacks_to_a_deterministic_canonical_wheel(installed, tmp_path):
    first, revision = downstream_runtime._materialize_installed_wheel(tmp_path / "cache-a")
    second, _ = downstream_runtime._materialize_installed_wheel(tmp_path / "cache-b")
    assert first.read_bytes() == second.read_bytes()
    digest = hashlib.sha256(first.read_bytes()).hexdigest()
    assert first.parent.name == digest and first.name == "nf_rna-1.2.0-py3-none-any.whl"
    assert revision is None

    with zipfile.ZipFile(first) as wheel:
        names = wheel.namelist()
        assert {info.date_time for info in wheel.infolist()} == {downstream_runtime.CANONICAL_WHEEL_TIMESTAMP}
        assert {info.compress_type for info in wheel.infolist()} == {zipfile.ZIP_STORED}
        record = wheel.read(f"{INFO}/RECORD").decode().splitlines()
        for line in record[:-1]:
            path, digest_field, size = line.split(",")
            assert digest_field == _record_hash(wheel.read(path)) and int(size) == len(wheel.read(path))
    assert names[-1] == f"{INFO}/RECORD"
    assert not any("__pycache__" in name or name.startswith("../") for name in names)
    assert not {f"{INFO}/INSTALLER", f"{INFO}/REQUESTED", f"{INFO}/direct_url.json"} & set(names)
    for bundled in ("rnaseq/workflow_support.py", "rnaseq/r/l2_analysis.R", "rnaseq/workflows/main.nf",
                    "rnaseq/workflows/hisat2_featurecounts.nf", "rnaseq/workflows/nextflow.config",
                    "rnaseq/workflows/envs/hisat2-featurecounts-linux-64.yml",
                    "rnaseq/workflows/envs/locks/nf-rna-downstream-linux-64.lock.yml",
                    "rnaseq/workflows/envs/locks/SHA256SUMS"):
        assert bundled in names
    assert (installed / "rnaseq" / "workflows" / "main.nf").read_bytes() == (ROOT / "workflow" / "main.nf").read_bytes()


def test_installer_metadata_does_not_change_runtime_identity(tmp_path, monkeypatch):
    digests = []
    for name, direct_url in (("a", {"url": "file:///x.whl", "archive_info": {}}), ("b", {"url": "https://pypi.org/x.whl", "archive_info": {}})):
        distribution = _install(tmp_path / name, direct_url=direct_url)
        monkeypatch.setattr(downstream_runtime, "_installed_distribution", lambda d=distribution: d)
        wheel, _ = downstream_runtime._materialize_installed_wheel(tmp_path / f"cache-{name}")
        digests.append(hashlib.sha256(wheel.read_bytes()).hexdigest())
    assert digests[0] == digests[1]


def test_vcs_install_records_its_commit_as_source_revision(tmp_path, monkeypatch):
    distribution = _install(tmp_path / "site", direct_url={"url": "https://github.com/x/nf-rna.git", "vcs_info": {"vcs": "git", "commit_id": "157410a4"}})
    monkeypatch.setattr(downstream_runtime, "_installed_distribution", lambda: distribution)
    assert downstream_runtime._materialize_installed_wheel(tmp_path / "cache")[1] == "157410a4"


def test_locally_modified_installation_is_refused(installed, tmp_path):
    script = installed / "rnaseq" / "r" / "l2_analysis.R"
    script.write_text(script.read_text(encoding="utf-8") + "# edited in site-packages\n", encoding="utf-8")
    with pytest.raises(ExecutionPreflightError, match="modified after installation: rnaseq/r/l2_analysis.R"):
        downstream_runtime._materialize_installed_wheel(tmp_path / "cache")
    ready, detail = downstream_runtime.wheel_source_status()
    assert not ready and "modified after installation" in detail


def test_incomplete_installation_is_refused(installed, tmp_path):
    (installed / "rnaseq" / "workflows" / "main.nf").unlink()
    with pytest.raises(ExecutionPreflightError, match="Installed nf-rna file is missing: rnaseq/workflows/main.nf"):
        downstream_runtime._materialize_installed_wheel(tmp_path / "cache")


def test_installed_nf_rna_provisions_and_reuses_the_locked_prefix(installed, monkeypatch, tmp_path):
    monkeypatch.setenv("RNASEQ_RUNTIME_ROOT", str(tmp_path / "runtime"))
    conda = FakeConda(monkeypatch)
    runtime = downstream_runtime.ensure_downstream_runtime()
    identity = runtime.identity()
    assert identity["wheel"]["origin"] == downstream_runtime.WHEEL_ORIGIN_INSTALLED
    assert identity["lock"]["sha256"] == "18197d49b7087d39c3a8e3aef8a16d7bd65af7e56e07d87aa47be656b219617d"
    assert runtime.prefix.name == "nf-rna-downstream-linux-64-18197d49b7087d39"
    marker = json.loads((runtime.prefix / "runtime" / "nf-rna-wheel.json").read_text(encoding="utf-8"))
    assert marker == {"filename": runtime.wheel_filename, "sha256": runtime.wheel_sha256, "origin": "installed-distribution"}
    # The R scripts inventoried inside the prefix are exactly the bundled source scripts.
    assert {item["name"]: item["sha256"] for item in runtime.r_scripts} == {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in (ROOT / "src" / "rnaseq" / "r").glob("*.R")
    }
    assert conda.created() == 1

    again = downstream_runtime.ensure_downstream_runtime()
    assert conda.created() == 1
    assert again.identity() == identity


def test_changed_release_is_rejected_by_the_cached_prefix(installed, monkeypatch, tmp_path):
    monkeypatch.setenv("RNASEQ_RUNTIME_ROOT", str(tmp_path / "runtime"))
    conda = FakeConda(monkeypatch)
    runtime = downstream_runtime.ensure_downstream_runtime()
    marker = runtime.prefix / "runtime" / "nf-rna-wheel.json"
    recorded = marker.read_bytes()

    # A different (properly installed) release: new content with a consistent RECORD.
    support = installed / "rnaseq" / "workflow_support.py"
    support.write_text(support.read_text(encoding="utf-8") + "\n# next release\n", encoding="utf-8")
    _write_record(installed)
    with pytest.raises(ExecutionPreflightError, match="wheel identity does not match") as error:
        downstream_runtime.ensure_downstream_runtime()
    assert runtime.wheel_sha256 in str(error.value)
    assert marker.read_bytes() == recorded
    assert conda.created() == 1


def test_prefix_without_identity_marker_is_never_reused(installed, monkeypatch, tmp_path):
    monkeypatch.setenv("RNASEQ_RUNTIME_ROOT", str(tmp_path / "runtime"))
    FakeConda(monkeypatch)
    runtime = downstream_runtime.ensure_downstream_runtime()
    (runtime.prefix / "runtime" / "nf-rna-wheel.json").unlink()
    with pytest.raises(ExecutionPreflightError, match="lacks its wheel identity marker"):
        downstream_runtime.ensure_downstream_runtime()


def test_source_checkout_development_mode_is_retained(monkeypatch, tmp_path):
    assert downstream_runtime._source_checkout() == ROOT.resolve()
    built = tmp_path / "nf_rna-1.2.0-py3-none-any.whl"
    monkeypatch.setattr(downstream_runtime, "_build_wheel", lambda _cache: (built, "abc+dirty"))
    monkeypatch.setattr(downstream_runtime, "_materialize_installed_wheel", lambda _cache: pytest.fail("source mode re-packed"))
    assert downstream_runtime._runtime_wheel(tmp_path) == (built, "abc+dirty", downstream_runtime.WHEEL_ORIGIN_SOURCE)


def test_non_editable_install_is_detected_in_a_real_interpreter(tmp_path):
    """A fresh interpreter importing nf-rna from site-packages needs no checkout or build tool."""

    _install(tmp_path / "site")
    probe = (
        "import json, pathlib, rnaseq\n"
        "from rnaseq import downstream_runtime as d\n"
        "print(json.dumps({'package': str(pathlib.Path(rnaseq.__file__).parent), 'checkout': str(d._source_checkout()),\n"
        "                  'status': d.wheel_source_status(), 'owner': str(d._installed_distribution()._path)}))\n"
    )
    environment = {**os.environ, "PYTHONPATH": str(tmp_path / "site")}
    result = subprocess.run([sys.executable, "-c", probe], cwd=tmp_path, env=environment, capture_output=True, text=True, check=True)
    observed = json.loads(result.stdout)
    assert observed["package"] == str((tmp_path / "site" / "rnaseq").resolve())
    assert observed["checkout"] == "None"
    assert observed["status"][0] is True and "installed distribution nf_rna-1.2.0" in observed["status"][1]
    # Another nf-rna installed elsewhere on sys.path must not be mistaken for the running one.
    assert observed["owner"] == str(tmp_path / "site" / INFO)
