"""Every run must be attributable to one exact, committed nf-rna source revision."""

from __future__ import annotations

import json

import pytest

from rnaseq import downstream_runtime
from rnaseq.errors import ExecutionPreflightError
from rnaseq.execution import downstream_runtime_checks


pytestmark = pytest.mark.real_source_revision
COMMIT = "62e2ce4e02989a887455bfa9eecfa1ba3110c99d"


class _Distribution:
    def __init__(self, direct_url: dict | None):
        self._direct_url = direct_url

    def read_text(self, name: str) -> str | None:
        assert name == "direct_url.json"
        return None if self._direct_url is None else json.dumps(self._direct_url)


def _installed(monkeypatch, direct_url: dict | None) -> None:
    monkeypatch.setattr(downstream_runtime, "_source_checkout", lambda: None)
    monkeypatch.setattr(downstream_runtime, "_installed_distribution", lambda: _Distribution(direct_url))


def test_git_install_is_identified_by_the_commit_pip_recorded(monkeypatch):
    _installed(monkeypatch, {"url": "https://github.com/2002brian/nf-rna.git", "vcs_info": {"vcs": "git", "commit_id": COMMIT, "requested_revision": "v1.3.0"}})
    assert downstream_runtime.runtime_source_revision() == COMMIT
    assert downstream_runtime.require_identified_source() == COMMIT


@pytest.mark.parametrize(
    "direct_url",
    [None, {"url": "file:///home/user/nf-rna", "dir_info": {}}, {"url": "file:///tmp/nf_rna-1.3.0-py3-none-any.whl", "archive_info": {}}],
    ids=["no-direct-url", "directory-install", "wheel-install"],
)
def test_install_without_a_recorded_commit_is_refused(monkeypatch, direct_url):
    _installed(monkeypatch, direct_url)
    assert downstream_runtime.runtime_source_revision() is None
    with pytest.raises(ExecutionPreflightError, match="source revision is unavailable.*git\\+https"):
        downstream_runtime.require_identified_source()


def test_source_checkout_reports_head_and_refuses_uncommitted_changes(monkeypatch, tmp_path):
    monkeypatch.setattr(downstream_runtime, "_source_checkout", lambda: tmp_path)
    monkeypatch.setattr(downstream_runtime, "_source_revision", lambda _root: COMMIT)
    assert downstream_runtime.require_identified_source() == COMMIT
    monkeypatch.setattr(downstream_runtime, "_source_revision", lambda _root: COMMIT + "+dirty")
    with pytest.raises(ExecutionPreflightError, match="uncommitted changes"):
        downstream_runtime.require_identified_source()


def test_run_preflight_refuses_an_unidentified_source_before_any_environment_work(monkeypatch):
    _installed(monkeypatch, {"url": "file:///home/user/nf-rna", "dir_info": {}})
    monkeypatch.setattr(downstream_runtime, "_lock_checksum", lambda _lock: pytest.fail("lock checked before identity"))
    with pytest.raises(ExecutionPreflightError, match="source revision is unavailable"):
        downstream_runtime.downstream_runtime_preflight()


def test_doctor_reports_the_source_revision(monkeypatch, tmp_path):
    monkeypatch.setenv("RNASEQ_RUNTIME_ROOT", str(tmp_path / "runtime"))
    _installed(monkeypatch, {"url": "https://github.com/2002brian/nf-rna.git", "vcs_info": {"vcs": "git", "commit_id": COMMIT}})
    monkeypatch.setattr(downstream_runtime, "wheel_source_status", lambda: (True, "installed"))
    checks = {check.name: check for check in downstream_runtime_checks()}
    assert checks["nf-rna source revision"].state == "FOUND" and COMMIT in checks["nf-rna source revision"].detail
    _installed(monkeypatch, {"url": "file:///home/user/nf-rna", "dir_info": {}})
    checks = {check.name: check for check in downstream_runtime_checks()}
    assert checks["nf-rna source revision"].state == "NOT FOUND"


def test_prefix_is_unique_per_lock_and_wheel(tmp_path):
    first = downstream_runtime.runtime_prefix(tmp_path, "linux-64", "a" * 64, "b" * 64)
    assert first == tmp_path / "prefixes" / f"nf-rna-downstream-linux-64-{'a' * 16}-{'b' * 16}"
    assert downstream_runtime.runtime_prefix(tmp_path, "linux-64", "a" * 64, "c" * 64) != first
