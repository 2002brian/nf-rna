"""First-party Docker image identity: version match, release revision and frozen digest."""

from __future__ import annotations

import json
from datetime import datetime
from types import SimpleNamespace

import pytest
import yaml

from rnaseq import service
from rnaseq.errors import ExecutionPreflightError
from rnaseq.execution import RuntimeCheck, check_container_runtime
from rnaseq.models import DEFAULT_EXECUTION_IMAGE, PIPELINE_VERSION
from rnaseq.planner import generate_plan
from rnaseq.service import (
    ExecutionImage, create_case_run, execute_service_run, freeze_case_inputs, prepare_service_run,
    resolve_downstream_inputs, resolve_execution_image, write_downstream_runtime_config,
)
from rnaseq.validators import validate_project
from rnaseq.workflow_support import l1_config, l2_config


REPOSITORY = DEFAULT_EXECUTION_IMAGE.rsplit(":", 1)[0]
DIGEST = f"{REPOSITORY}@sha256:" + "a" * 64
RELEASE = "1" * 40


def _probe(monkeypatch, stdout: str) -> RuntimeCheck:
    monkeypatch.setattr("rnaseq.execution.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.execution._run_capture", lambda _arguments: SimpleNamespace(returncode=0, stdout=stdout, stderr=""))
    return check_container_runtime(DEFAULT_EXECUTION_IMAGE)


def _inspected(monkeypatch, *, digests: list[str], revision: str | None = RELEASE, image_id: str | None = "sha256:" + "a" * 64):
    labels = {"org.opencontainers.image.revision": revision} if revision is not None else {}
    monkeypatch.setattr(service, "inspect_container_image", lambda image: {
        "reference": image, "image_id": image_id, "repo_digests": digests, "labels": labels,
    })


# A. Version -----------------------------------------------------------------

def test_image_reporting_the_host_version_is_accepted(monkeypatch):
    result = _probe(monkeypatch, f"nf-rna-version={PIPELINE_VERSION}\n")
    assert result.state == "FOUND"
    assert f"nf-rna {PIPELINE_VERSION}" in result.detail


def test_older_image_version_is_incompatible(monkeypatch):
    result = _probe(monkeypatch, "nf-rna-version=1.1.2\n")
    assert result.state == "INCOMPATIBLE"
    assert "contains nf-rna 1.1.2" in result.detail and f"this CLI is nf-rna {PIPELINE_VERSION}" in result.detail


@pytest.mark.parametrize("stdout", ["", "nf-rna-version=\n", f"nf-rna-version={PIPELINE_VERSION}\nnf-rna-version=1.1.2\n"])
def test_missing_or_malformed_image_version_is_rejected(monkeypatch, stdout):
    result = _probe(monkeypatch, stdout)
    assert result.state == "NOT FOUND"
    assert "did not report its installed nf-rna version" in result.detail


# B. Revision ----------------------------------------------------------------

def test_production_accepts_clean_checkout_at_the_image_revision(monkeypatch):
    _inspected(monkeypatch, digests=[DIGEST])
    monkeypatch.setattr(service, "_host_source_checkout", lambda: ("/src", RELEASE, False))
    assert resolve_execution_image(DEFAULT_EXECUTION_IMAGE, production=True).reference == DIGEST


def test_production_accepts_installed_cli_without_a_checkout(monkeypatch):
    _inspected(monkeypatch, digests=[DIGEST])
    monkeypatch.setattr(service, "_host_source_checkout", lambda: (None, None, False))
    assert resolve_execution_image(DEFAULT_EXECUTION_IMAGE, production=True).source_revision == RELEASE


@pytest.mark.parametrize(("head", "dirty"), [("2" * 40, False), (RELEASE, True)])
def test_production_rejects_mismatched_or_dirty_host_checkout(monkeypatch, head, dirty):
    _inspected(monkeypatch, digests=[DIGEST])
    monkeypatch.setattr(service, "_host_source_checkout", lambda: ("/src", head, dirty))
    with pytest.raises(ExecutionPreflightError, match="requires a clean checkout at the execution image revision"):
        resolve_execution_image(DEFAULT_EXECUTION_IMAGE, production=True)


@pytest.mark.parametrize("revision", [None, "unlabeled-development-image", "abc123", RELEASE + "+dirty", "A" * 40])
def test_production_rejects_missing_malformed_or_dirty_image_revision(monkeypatch, revision):
    _inspected(monkeypatch, digests=[DIGEST], revision=revision)
    monkeypatch.setattr(service, "_host_source_checkout", lambda: (None, None, False))
    with pytest.raises(ExecutionPreflightError, match="requires a released execution image"):
        resolve_execution_image(DEFAULT_EXECUTION_IMAGE, production=True)


# C. Digest ------------------------------------------------------------------

def test_repo_digest_for_the_requested_repository_is_chosen_deterministically(monkeypatch):
    other = "example.org/other@sha256:" + "0" * 64
    second = f"{REPOSITORY}@sha256:" + "b" * 64
    _inspected(monkeypatch, digests=[second, other, DIGEST])
    monkeypatch.setattr(service, "_host_source_checkout", lambda: (None, None, False))
    assert resolve_execution_image(DEFAULT_EXECUTION_IMAGE, production=True).repo_digest == DIGEST
    # A requested digest reference is kept as-is.
    assert resolve_execution_image(second, production=True).reference == second


@pytest.mark.parametrize(("digests", "image_id"), [([], "sha256:" + "a" * 64), (["example.org/other@sha256:" + "0" * 64], "sha256:" + "a" * 64), ([DIGEST], None)])
def test_production_without_immutable_identity_is_rejected(monkeypatch, digests, image_id):
    _inspected(monkeypatch, digests=digests, image_id=image_id)
    with pytest.raises(ExecutionPreflightError, match="observed immutable execution image ID/digest"):
        resolve_execution_image(DEFAULT_EXECUTION_IMAGE, production=True)


def test_preflight_probes_and_freezes_the_resolved_digest(monkeypatch, project_factory, production_capable_execution_capacity):
    report = validate_project(project_factory())
    generate_plan(report)
    _inspected(monkeypatch, digests=[DIGEST], revision="unlabeled-container-image")
    probed: list[str] = []
    monkeypatch.setattr(service, "check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "test"))
    monkeypatch.setattr(service, "check_docker", lambda: RuntimeCheck("Docker", "FOUND", "test"))
    monkeypatch.setattr(service, "check_container_runtime", lambda image: probed.append(image) or RuntimeCheck("image", "FOUND", "test"))
    _resources, image = prepare_service_run(report, profile="local")
    assert probed == [DIGEST]
    assert (image.requested, image.reference, image.nf_rna_version) == (DEFAULT_EXECUTION_IMAGE, DIGEST, PIPELINE_VERSION)

    run = create_case_run(report, "CASE-IDENTITY", moment=datetime(2026, 9, 24, 12, 0, 0))
    frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"], execution_image=image)
    contract = json.loads(frozen.contract.read_text(encoding="utf-8"))["execution"]
    manifest = yaml.safe_load((run.run_dir / "frozen" / "execution_manifest.yaml").read_text(encoding="utf-8"))
    assert contract == {
        "image": DEFAULT_EXECUTION_IMAGE, "resolved_image": DIGEST, "image_id": "sha256:" + "a" * 64,
        "nf_rna_version": PIPELINE_VERSION, "source_revision": "unlabeled-container-image",
    }
    assert (manifest["execution_image"], manifest["resolved_execution_image"], manifest["execution_nf_rna_version"]) == (
        DEFAULT_EXECUTION_IMAGE, DIGEST, PIPELINE_VERSION,
    )
    runtime_config = write_downstream_runtime_config(run, image.reference)
    assert runtime_config.read_text(encoding="utf-8") == f'params.first_party_image = "{DIGEST}"\n'


def test_incompatible_image_fails_before_any_run_is_allocated(monkeypatch, project_factory, production_capable_execution_capacity):
    root = project_factory()
    report = validate_project(root)
    generate_plan(report)
    _inspected(monkeypatch, digests=[DIGEST])
    monkeypatch.setattr(service, "check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "test"))
    monkeypatch.setattr(service, "check_docker", lambda: RuntimeCheck("Docker", "FOUND", "test"))
    monkeypatch.setattr(service, "check_container_runtime", lambda _image: RuntimeCheck("image", "INCOMPATIBLE", "contains nf-rna 1.1.2"))
    monkeypatch.setattr(service, "_run_command", lambda *_args, **_kwargs: pytest.fail("Nextflow must not run"))
    with pytest.raises(ExecutionPreflightError, match="contains nf-rna 1.1.2"):
        execute_service_run(report, case_id="CASE-OLD-IMAGE")
    assert not (root / "runs").exists()


def test_retry_verifies_the_source_runs_frozen_digest(monkeypatch, project_factory, production_capable_execution_capacity):
    report = validate_project(project_factory())
    _inspected(monkeypatch, digests=[DIGEST])
    probed: list[str] = []
    monkeypatch.setattr(service, "check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "test"))
    monkeypatch.setattr(service, "check_docker", lambda: RuntimeCheck("Docker", "FOUND", "test"))
    monkeypatch.setattr(service, "check_container_runtime", lambda image: probed.append(image) or RuntimeCheck("image", "FOUND", "test"))
    _resources, image = service._prepare_retry_runtime(report, "local", DIGEST)
    assert probed == [DIGEST] and image.reference == DIGEST


# D. Development -------------------------------------------------------------

def test_non_production_image_without_repo_digest_runs_by_its_requested_reference(monkeypatch):
    _inspected(monkeypatch, digests=[], revision="unlabeled-development-image")
    image = resolve_execution_image("nf-rna:dev", production=False)
    assert (image.reference, image.repo_digest, image.source_revision) == ("nf-rna:dev", None, "unlabeled-development-image")


# E. Scientific configuration ------------------------------------------------

def test_image_identity_does_not_change_scientific_configuration(monkeypatch, project_factory, tmp_path):
    report = validate_project(project_factory())
    generate_plan(report)
    _inspected(monkeypatch, digests=[])
    images = (
        ExecutionImage(DEFAULT_EXECUTION_IMAGE, DEFAULT_EXECUTION_IMAGE, None, None, "unlabeled-container-image"),
        ExecutionImage(DEFAULT_EXECUTION_IMAGE, DIGEST, DIGEST, "sha256:" + "a" * 64, RELEASE, PIPELINE_VERSION),
    )
    views = []
    for index, image in enumerate(images):
        run = create_case_run(report, f"CASE-SCIENCE-{index}", moment=datetime(2026, 9, 24, 12, 0, 0))
        frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"], execution_image=image)
        contract = json.loads(frozen.contract.read_text(encoding="utf-8"))
        inputs = resolve_downstream_inputs(run).root
        configs = [l1_config(frozen.contract, inputs, tmp_path / "l1"), l2_config(frozen.contract, inputs, tmp_path / "l1", tmp_path / "l2")]
        # Per-run paths and runtime provenance differ by design; everything scientific must not.
        for config in configs:
            for key in ("metadata", "counts", "runtime"):
                config.pop(key)
        contract["source"].pop("counts")
        views.append(({key: contract[key] for key in ("source", "analysis_level", "analysis", "design", "annotation")}, configs))
    assert views[0] == views[1]
