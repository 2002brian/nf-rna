"""Audit D-07 regression: the delivered run_state is the final run state.

``assemble_delivery`` used to copy ``run_state.json`` into
``delivery/methods_and_versions`` while the run was still ``RUNNING`` (phase
``delivery``), so the client package recorded a non-final state with
``completed_at: null``.  The delivered state is now the SUCCESS state that is
recorded once the package verifies, with the same ``completed_at``.
"""

from __future__ import annotations

import json

import pytest

from rnaseq.planner import generate_plan
from rnaseq.service import execute_service_run, verify_delivery_manifest
from rnaseq.validators import validate_project
from test_runtime_dispatch import LINUX, _fake_nextflow, _guard_conda, _host, mocked_downstream_runtime  # noqa: F401


pytestmark = pytest.mark.usefixtures("production_capable_execution_capacity")


def _successful_raw_run(monkeypatch, tmp_path, project_factory):
    report = validate_project(project_factory())
    generate_plan(report)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    _host(monkeypatch, *LINUX)
    _guard_conda(monkeypatch)
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow([]))
    return execute_service_run(report, case_id="CASE-AUDIT")


def _delivered_state(run):
    (path,) = (run.run_dir / "delivery" / "methods_and_versions").glob("run_state_*.json")
    return json.loads(path.read_text(encoding="utf-8"))


def test_run_is_success_and_delivery_manifest_verifies(monkeypatch, tmp_path, project_factory, mocked_downstream_runtime):
    run = _successful_raw_run(monkeypatch, tmp_path, project_factory)
    assert json.loads(run.state_path.read_text(encoding="utf-8"))["status"] == "SUCCESS"
    assert verify_delivery_manifest(run.run_dir / "delivery") == ()


def test_delivered_run_state_is_final(monkeypatch, tmp_path, project_factory, mocked_downstream_runtime):
    run = _successful_raw_run(monkeypatch, tmp_path, project_factory)
    delivered = _delivered_state(run)
    assert delivered["status"] == "SUCCESS"
    assert delivered["completed_at"] is not None


def test_delivered_run_state_is_exactly_the_final_run_state(monkeypatch, tmp_path, project_factory, mocked_downstream_runtime):
    run = _successful_raw_run(monkeypatch, tmp_path, project_factory)
    final = json.loads(run.state_path.read_text(encoding="utf-8"))
    assert _delivered_state(run) == final
    assert final["phase"] == "delivery" and final["delivery"] == str(run.run_dir / "delivery")

