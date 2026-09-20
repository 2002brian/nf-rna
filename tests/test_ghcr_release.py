"""Behavioral safeguards for release-only GHCR publication."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from urllib.error import URLError

import pytest


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ghcr_release  # noqa: E402


def _error(code: str) -> bytes:
    return ('{"errors":[{"code":"' + code + '"}]}').encode()


@pytest.mark.parametrize(
    ("tag", "package", "container", "stable"),
    [
        ("v1.1.0", "1.1.0", "1.1.0", True),
        ("v1.1.1", "1.1.1", "1.1.1", True),
        ("v2.0.0", "2.0.0", "2.0.0", True),
        ("v1.1.0-rc1", "1.1.0rc1", "1.1.0-rc1", False),
        ("v1.1.0-a2", "1.1.0a2", "1.1.0-a2", False),
        ("v1.1.0-b3", "1.1.0b3", "1.1.0-b3", False),
    ],
)
def test_canonical_release_tags_have_one_package_and_container_identity(tag, package, container, stable):
    identity = ghcr_release.parse_release_tag(tag)

    assert identity.package_version == package
    assert identity.container_version == container
    assert identity.stable is stable


@pytest.mark.parametrize("tag", ["v1.1.0rc1", "v1.1.0-rc0", "1.1.0", "v1.1", "v1.1.0-dev1"])
def test_ambiguous_or_malformed_release_tags_are_rejected(tag):
    with pytest.raises(ValueError):
        ghcr_release.parse_release_tag(tag)


def test_release_identity_rejects_package_or_checkout_mismatch():
    with pytest.raises(ValueError, match="Package version"):
        ghcr_release.verify_release_identity("v1.1.0-rc1", package_version="1.1.0")
    with pytest.raises(ValueError, match="Checkout"):
        ghcr_release.verify_release_identity(
            "v1.1.0", package_version="1.1.0", tag_revision="a" * 40, checkout_revision="b" * 40
        )


@pytest.mark.parametrize("kind", ["annotated", "lightweight"])
def test_tag_resolution_uses_commit_dereference_for_annotated_and_lightweight_tags(kind):
    seen: list[list[str]] = []

    def fake_run(command, **kwargs):
        seen.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="a" * 40 + "\n", stderr="")

    assert ghcr_release.resolve_tag_commit("v1.1.0", fake_run) == "a" * 40
    assert seen == [["git", "rev-parse", "--verify", "refs/tags/v1.1.0^{commit}"]]


@pytest.mark.parametrize(
    ("status", "body", "decision", "reason"),
    [
        (200, b"{}", ghcr_release.RegistryDecision.PRESENT, "manifest_present"),
        (404, _error("MANIFEST_UNKNOWN"), ghcr_release.RegistryDecision.ABSENT, "manifest_absent"),
        (404, _error("NAME_UNKNOWN"), ghcr_release.RegistryDecision.ABSENT, "package_absent_first_publication"),
        (401, _error("UNAUTHORIZED"), ghcr_release.RegistryDecision.ERROR, "unexpected_http_401"),
        (403, _error("DENIED"), ghcr_release.RegistryDecision.ERROR, "unexpected_http_403"),
        (429, _error("TOOMANYREQUESTS"), ghcr_release.RegistryDecision.ERROR, "unexpected_http_429"),
        (500, _error("UNKNOWN"), ghcr_release.RegistryDecision.ERROR, "unexpected_http_500"),
        (418, b"{}", ghcr_release.RegistryDecision.ERROR, "unexpected_http_418"),
        (404, b"not json", ghcr_release.RegistryDecision.ERROR, "malformed_404_response"),
        (404, _error("BLOB_UNKNOWN"), ghcr_release.RegistryDecision.ERROR, "unexpected_404_response"),
    ],
)
def test_registry_decision_table_is_structured_and_fail_closed(status, body, decision, reason):
    result = ghcr_release.assess_manifest_response(status, body)

    assert result.decision is decision
    assert result.reason == reason


class _Response:
    def __init__(self, status: int, headers: dict[str, str], body: bytes):
        self.status = status
        self.headers = headers
        self._body = body

    def read(self) -> bytes:
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def test_registry_check_follows_bearer_challenge_and_sends_manifest_accept_header():
    calls = []
    responses = iter(
        [
            _Response(401, {"WWW-Authenticate": 'Bearer realm="https://ghcr.io/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"'}, b""),
            _Response(200, {}, b'{"token":"registry-token"}'),
            _Response(404, {}, _error("MANIFEST_UNKNOWN")),
        ]
    )

    def opener(request, timeout):
        calls.append(request)
        return next(responses)

    result = ghcr_release.check_registry_tag("ghcr.io", "2002brian/nf-rna", "1.1.0", "actor", "token", opener)

    assert result.decision is ghcr_release.RegistryDecision.ABSENT
    assert calls[0].get_header("Accept") == ghcr_release.OCI_MANIFEST_ACCEPT
    assert calls[0].get_header("Authorization", "").startswith("Basic ")
    assert calls[2].get_header("Authorization") == "Bearer registry-token"


@pytest.mark.parametrize("failure", [TimeoutError(), URLError("network unavailable")])
def test_registry_transport_failures_cannot_permit_publication(failure):
    def opener(request, timeout):
        raise failure

    result = ghcr_release.check_registry_tag("ghcr.io", "2002brian/nf-rna", "1.1.0", "actor", "token", opener)

    assert result.decision is ghcr_release.RegistryDecision.ERROR
    assert result.reason.startswith("registry_communication_error:")


def test_image_validation_requires_configured_and_effective_non_root_and_all_assets(monkeypatch):
    calls = []

    def fake_docker(*args, input_text=None):
        calls.append(args)
        if args[:4] == ("image", "inspect", "--format", "{{.Config.User}}"):
            return "mambauser"
        if args[:3] == ("image", "inspect", "--format"):
            if "source" in args[3]:
                return "https://github.com/2002brian/nf-rna"
            if "revision" in args[3]:
                return "a" * 40
            return "1.1.0"
        if args[-2:] == ("id", "-u"):
            return "1000"
        if args[-2:] == ("rnaseq", "--version"):
            return "1.1.0"
        if args[-1].startswith("from rnaseq import"):
            return "1.1.0"
        return ""

    monkeypatch.setattr(ghcr_release, "_docker", fake_docker)
    ghcr_release.validate_image("image", "1.1.0", "1.1.0", "a" * 40)

    checked_r_assets = {args[-1] for args in calls if args[:4] == ("run", "--rm", "image", "test")}
    assert checked_r_assets == {f"/opt/nf-rna/r/{asset}" for asset in ghcr_release.REQUIRED_R_ASSETS}


def test_image_validation_rejects_effective_root(monkeypatch):
    def fake_docker(*args, input_text=None):
        if args[:4] == ("image", "inspect", "--format", "{{.Config.User}}"):
            return "mambauser"
        if args[-2:] == ("id", "-u"):
            return "0"
        raise AssertionError(f"Unexpected command after root check: {args}")

    monkeypatch.setattr(ghcr_release, "_docker", fake_docker)
    with pytest.raises(ValueError, match="fell back to root"):
        ghcr_release.validate_image("image", "1.1.0", "1.1.0", "a" * 40)


VALID_CHALLENGE = (
    'Bearer realm="https://ghcr.io/token",service="ghcr.io",'
    'scope="repository:2002brian/nf-rna:pull"'
)


@pytest.mark.parametrize(
    "challenge",
    [
        'Bearer realm="https://token.example/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io.example/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://subdomain.ghcr.io/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="http://ghcr.io/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io:444/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://user@ghcr.io/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io/unexpected",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io/token",service="other",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io/token",service="ghcr.io",scope="repository:other/image:pull"',
        'Bearer realm="https://ghcr.io/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull,push"',
        'Bearer service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io/token",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io/token",service="ghcr.io"',
        'Bearer realm="https://ghcr.io/token",realm="https://ghcr.io/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"',
        'Bearer realm="https://ghcr.io/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull" trailing',
        "Basic realm=ghcr.io",
    ],
)
def test_untrusted_or_malformed_bearer_challenges_are_rejected_before_token_request(challenge):
    with pytest.raises(ValueError):
        ghcr_release._trusted_token_url(challenge)


def test_exact_ghcr_bearer_challenge_constructs_only_the_trusted_token_endpoint():
    assert ghcr_release._trusted_token_url(VALID_CHALLENGE) == (
        "https://ghcr.io/token?service=ghcr.io&scope=repository%3A2002brian%2Fnf-rna%3Apull"
    )


def test_invalid_challenge_never_sends_basic_credentials_to_the_challenge_destination():
    calls = []

    def opener(request, timeout):
        calls.append(request)
        return _Response(
            401,
            {"WWW-Authenticate": 'Bearer realm="https://token.example/token",service="ghcr.io",scope="repository:2002brian/nf-rna:pull"'},
            b"",
        )

    result = ghcr_release.check_registry_tag("ghcr.io", "2002brian/nf-rna", "1.1.0", "actor", "token", opener)

    assert result.decision is ghcr_release.RegistryDecision.ERROR
    assert [request.full_url for request in calls] == ["https://ghcr.io/v2/2002brian/nf-rna/manifests/1.1.0"]
    assert calls[0].get_header("Authorization", "").startswith("Basic ")


@pytest.mark.parametrize(
    "location",
    [
        "https://other.example/manifest",
        "http://ghcr.io/v2/2002brian/nf-rna/manifests/1.1.0",
        "https://ghcr.io/unexpected",
        "https://ghcr.io/v2/2002brian/nf-rna/manifests/1.1.0",
    ],
)
def test_manifest_redirects_fail_closed_without_following_location(location):
    calls = []

    def opener(request, timeout):
        calls.append(request)
        return _Response(302, {"Location": location}, b"")

    result = ghcr_release.check_registry_tag("ghcr.io", "2002brian/nf-rna", "1.1.0", "actor", "token", opener)

    assert result.decision is ghcr_release.RegistryDecision.ERROR
    assert len(calls) == 1


def test_token_endpoint_redirect_fails_closed_without_a_bearer_manifest_request():
    calls = []
    responses = iter(
        [
            _Response(401, {"WWW-Authenticate": VALID_CHALLENGE}, b""),
            _Response(302, {"Location": "https://other.example/token"}, b""),
        ]
    )

    def opener(request, timeout):
        calls.append(request)
        return next(responses)

    result = ghcr_release.check_registry_tag("ghcr.io", "2002brian/nf-rna", "1.1.0", "actor", "token", opener)

    assert result.decision is ghcr_release.RegistryDecision.ERROR
    assert [request.full_url for request in calls] == [
        "https://ghcr.io/v2/2002brian/nf-rna/manifests/1.1.0",
        "https://ghcr.io/token?service=ghcr.io&scope=repository%3A2002brian%2Fnf-rna%3Apull",
    ]


def _image_docker(overrides: dict[str, object]):
    def fail(args):
        raise subprocess.CalledProcessError(1, ["docker", *args])

    def fake(*args, input_text=None):
        if args[:4] == ("image", "inspect", "--format", "{{.Config.User}}"):
            return str(overrides.get("configured_user", "mambauser"))
        if args[:3] == ("image", "inspect", "--format"):
            if "source" in args[3]:
                return str(overrides.get("source", "https://github.com/2002brian/nf-rna"))
            if "revision" in args[3]:
                return str(overrides.get("revision", "a" * 40))
            return str(overrides.get("version", "1.1.0"))
        if args[-2:] == ("id", "-u"):
            return str(overrides.get("effective_uid", "1000"))
        if args[-2:] == ("rnaseq", "--version"):
            return str(overrides.get("cli_version", "1.1.0"))
        if args[-1].startswith("from rnaseq import"):
            return str(overrides.get("package_version", "1.1.0"))
        if args[-2:] == ("report", "--help") and overrides.get("report_failure"):
            fail(args)
        if "required_workflow_assets" in args[-1] and overrides.get("workflow_asset_failure"):
            fail(args)
        if args[:4] == ("run", "--rm", "image", "test") and args[-1] in overrides.get("missing_r_assets", set()):
            fail(args)
        if args[:4] == ("run", "--rm", "image", "sh") and overrides.get("runtime_command_failure"):
            fail(args)
        if args[-1] == "import rnaseq.models, rnaseq.workflow_support" and overrides.get("python_dependency_failure"):
            fail(args)
        if args[:4] == ("run", "--rm", "image", "Rscript") and overrides.get("r_dependency_failure"):
            fail(args)
        return ""

    return fake


@pytest.mark.parametrize(
    ("label", "value"),
    [
        ("source", "https://wrong.example/nf-rna"),
        ("source", ""),
        ("revision", "b" * 40),
        ("revision", ""),
        ("revision", "unknown"),
        ("revision", "dirty"),
        ("revision", "development"),
        ("version", "1.1.1"),
        ("version", ""),
    ],
)
def test_image_validation_rejects_missing_or_mismatched_oci_labels(monkeypatch, label, value):
    monkeypatch.setattr(ghcr_release, "_docker", _image_docker({label: value}))

    with pytest.raises(ValueError, match="OCI label"):
        ghcr_release.validate_image("image", "1.1.0", "1.1.0", "a" * 40)


@pytest.mark.parametrize("configured_user", ["", "root", "0", "0:1000"])
def test_image_validation_rejects_empty_or_root_configured_user(monkeypatch, configured_user):
    monkeypatch.setattr(ghcr_release, "_docker", _image_docker({"configured_user": configured_user}))

    with pytest.raises(ValueError, match="configured to run as root"):
        ghcr_release.validate_image("image", "1.1.0", "1.1.0", "a" * 40)


@pytest.mark.parametrize("asset", ["main.nf", "nextflow.config"])
def test_image_validation_rejects_missing_canonical_workflow_assets(monkeypatch, asset):
    monkeypatch.setattr(ghcr_release, "_docker", _image_docker({"workflow_asset_failure": asset}))

    with pytest.raises(subprocess.CalledProcessError):
        ghcr_release.validate_image("image", "1.1.0", "1.1.0", "a" * 40)


@pytest.mark.parametrize("asset", ghcr_release.REQUIRED_R_ASSETS)
def test_image_validation_rejects_every_missing_required_r_asset(monkeypatch, asset):
    monkeypatch.setattr(ghcr_release, "_docker", _image_docker({"missing_r_assets": {f"/opt/nf-rna/r/{asset}"}}))

    with pytest.raises(subprocess.CalledProcessError):
        ghcr_release.validate_image("image", "1.1.0", "1.1.0", "a" * 40)


@pytest.mark.parametrize(
    "overrides",
    [
        {"report_failure": True},
        {"runtime_command_failure": True},
        {"python_dependency_failure": True},
        {"r_dependency_failure": True},
        {"cli_version": "1.1.1"},
        {"package_version": "1.1.1"},
    ],
)
def test_image_validation_rejects_report_dependency_and_version_failures(monkeypatch, overrides):
    monkeypatch.setattr(ghcr_release, "_docker", _image_docker(overrides))

    with pytest.raises((ValueError, subprocess.CalledProcessError)):
        ghcr_release.validate_image("image", "1.1.0", "1.1.0", "a" * 40)
