#!/usr/bin/env python3
"""Release-only GHCR publication safeguards for the nf-rna CI workflow.

This is deliberately CI tooling, not part of the installed scientific runtime.
It centralizes release identity, OCI registry decisions, and local-image checks
so the workflow and its regression tests exercise the same behavior.
"""

from __future__ import annotations

import argparse
import ast
import base64
import json
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Callable, Mapping
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener


RELEASE_TAG_PATTERN = re.compile(
    r"^v(?P<base>[0-9]+\.[0-9]+\.[0-9]+)(?:-(?P<stage>a|b|rc)(?P<serial>[1-9][0-9]*))?$"
)
OCI_MANIFEST_ACCEPT = ", ".join(
    (
        "application/vnd.oci.image.index.v1+json",
        "application/vnd.oci.image.manifest.v1+json",
        "application/vnd.docker.distribution.manifest.list.v2+json",
        "application/vnd.docker.distribution.manifest.v2+json",
    )
)
GHCR_REGISTRY = "ghcr.io"
GHCR_REPOSITORY = "2002brian/nf-rna"
GHCR_PULL_SCOPE = f"repository:{GHCR_REPOSITORY}:pull"
REQUIRED_R_ASSETS = (
    "l1_analysis.R",
    "l2_analysis.R",
    "go_analysis.R",
    "gsea_analysis.R",
    "kegg_analysis.R",
    "annotation_mapping_qc.R",
    "gsea_core_members.R",
    "gsea_term_filtering.R",
    "kegg_core_members.R",
    "ora_helpers.R",
    "provenance.R",
)
INVALID_REVISION_VALUES = {"", "unknown", "dirty", "development", "unlabeled-development-image"}


@dataclass(frozen=True)
class ReleaseIdentity:
    git_tag: str
    package_version: str
    container_version: str
    stable: bool


class RegistryDecision(str, Enum):
    ABSENT = "absent"
    PRESENT = "present"
    ERROR = "error"


@dataclass(frozen=True)
class RegistryResult:
    decision: RegistryDecision
    reason: str
    status: int | None = None


@dataclass(frozen=True)
class HttpResult:
    status: int
    headers: Mapping[str, str]
    body: bytes


class _RejectRedirects(HTTPRedirectHandler):
    """Turn every redirect into a response the registry policy can reject."""

    def redirect_request(self, request, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
        return None


_NO_REDIRECT_OPENER = build_opener(_RejectRedirects())


def _open_without_redirects(request: Request, timeout: int = 20):
    return _NO_REDIRECT_OPENER.open(request, timeout=timeout)


def parse_release_tag(tag: str) -> ReleaseIdentity:
    """Parse the one supported Git tag spelling into package and image forms."""

    match = RELEASE_TAG_PATTERN.fullmatch(tag)
    if not match:
        raise ValueError(
            "Release tag must be vX.Y.Z or vX.Y.Z-(a|b|rc)N; alternate PEP 440 "
            "spellings are deliberately rejected."
        )
    base = match.group("base")
    stage = match.group("stage")
    serial = match.group("serial")
    suffix = "" if stage is None else f"-{stage}{serial}"
    return ReleaseIdentity(
        git_tag=tag,
        package_version=base if stage is None else f"{base}{stage}{serial}",
        container_version=f"{base}{suffix}",
        stable=stage is None,
    )


def read_package_version(path: Path) -> str:
    """Read the literal single-source package version without importing project code."""

    module = ast.parse(path.read_text(encoding="utf-8"))
    values = [
        ast.literal_eval(node.value)
        for node in module.body
        if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name) and target.id == "__version__"
    ]
    if len(values) != 1 or not isinstance(values[0], str):
        raise ValueError(f"Could not read one literal __version__ from {path}")
    return values[0]


def resolve_tag_commit(tag: str, run: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run) -> str:
    """Dereference either an annotated or lightweight tag to its commit."""

    completed = run(
        ["git", "rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}"],
        check=True,
        text=True,
        capture_output=True,
    )
    revision = completed.stdout.strip()
    if not revision:
        raise ValueError(f"Tag {tag} did not resolve to a commit.")
    return revision


def verify_release_identity(
    tag: str,
    package_version: str | None = None,
    tag_revision: str | None = None,
    checkout_revision: str | None = None,
) -> ReleaseIdentity:
    identity = parse_release_tag(tag)
    if package_version is not None and package_version != identity.package_version:
        raise ValueError(
            f"Package version {package_version} does not match release tag {tag} "
            f"(expected {identity.package_version})."
        )
    if (tag_revision is None) != (checkout_revision is None):
        raise ValueError("Both tag and checkout revisions are required when verifying checkout identity.")
    if tag_revision is not None and tag_revision != checkout_revision:
        raise ValueError(f"Checkout ({checkout_revision}) does not match {tag} ({tag_revision}).")
    return identity


def assess_manifest_response(status: int, body: bytes) -> RegistryResult:
    """Classify an OCI manifest response solely by status and structured JSON."""

    if status == 200:
        return RegistryResult(RegistryDecision.PRESENT, "manifest_present", status)
    if status != 404:
        return RegistryResult(RegistryDecision.ERROR, f"unexpected_http_{status}", status)
    try:
        payload = json.loads(body.decode("utf-8"))
        codes = {error["code"] for error in payload["errors"] if isinstance(error, dict)}
    except (UnicodeDecodeError, json.JSONDecodeError, KeyError, TypeError):
        return RegistryResult(RegistryDecision.ERROR, "malformed_404_response", status)
    if codes == {"MANIFEST_UNKNOWN"}:
        return RegistryResult(RegistryDecision.ABSENT, "manifest_absent", status)
    if codes == {"NAME_UNKNOWN"}:
        # GHCR may not create a package until its first manifest is pushed.
        return RegistryResult(RegistryDecision.ABSENT, "package_absent_first_publication", status)
    return RegistryResult(RegistryDecision.ERROR, "unexpected_404_response", status)


def _basic_authorization(username: str, token: str) -> str:
    encoded = base64.b64encode(f"{username}:{token}".encode("utf-8")).decode("ascii")
    return f"Basic {encoded}"


def _open(request: Request, opener: Callable[..., object]) -> HttpResult:
    def normalized_headers(headers: object) -> dict[str, str]:
        """Store HTTP field names canonically; HTTP header names are case-insensitive."""

        return {str(name).lower(): str(value) for name, value in headers.items()}  # type: ignore[union-attr]

    try:
        with opener(request, timeout=20) as response:
            return HttpResult(response.status, normalized_headers(response.headers), response.read())
    except HTTPError as error:
        return HttpResult(error.code, normalized_headers(error.headers), error.read())


def _parse_bearer_challenge(challenge: str) -> dict[str, str]:
    """Parse a strict, quoted OCI Bearer challenge without duplicate fields."""

    match = re.fullmatch(r"Bearer[ \t]+(.+)", challenge, flags=re.IGNORECASE)
    if not match:
        raise ValueError("Registry did not provide a supported Bearer authentication challenge.")
    remainder = match.group(1)
    position = 0
    attributes: dict[str, str] = {}
    while position < len(remainder):
        field = re.match(r"[A-Za-z]+", remainder[position:])
        if not field:
            raise ValueError("Registry Bearer challenge has malformed parameter syntax.")
        name = field.group(0).lower()
        position += len(field.group(0))
        if position >= len(remainder) or remainder[position] != "=":
            raise ValueError("Registry Bearer challenge has malformed parameter syntax.")
        position += 1
        if position >= len(remainder) or remainder[position] != '"':
            raise ValueError("Registry Bearer challenge parameters must be quoted.")
        position += 1
        end = remainder.find('"', position)
        if end < 0 or "\\" in remainder[position:end]:
            raise ValueError("Registry Bearer challenge has malformed quoting or escaping.")
        value = remainder[position:end]
        if name in attributes or name not in {"realm", "service", "scope"}:
            raise ValueError("Registry Bearer challenge has duplicate or unsupported parameters.")
        attributes[name] = value
        position = end + 1
        if position == len(remainder):
            break
        if remainder[position] != ",":
            raise ValueError("Registry Bearer challenge has malformed parameter separation.")
        position += 1
        if position >= len(remainder) or remainder[position].isspace():
            raise ValueError("Registry Bearer challenge has malformed parameter separation.")
    if set(attributes) != {"realm", "service", "scope"}:
        raise ValueError("Registry Bearer challenge must provide exactly realm, service, and scope.")
    return attributes


def _trusted_token_url(challenge: str) -> str:
    """Validate a GHCR challenge, then construct its token URL from trusted values."""

    attributes = _parse_bearer_challenge(challenge)
    realm = urlparse(attributes["realm"])
    try:
        realm_port = realm.port
    except ValueError as error:
        raise ValueError("Registry Bearer realm has an invalid port.") from error
    if (
        realm.scheme != "https"
        or realm.hostname != GHCR_REGISTRY
        or realm_port not in {None, 443}
        or realm.path != "/token"
        or realm.params
        or realm.query
        or realm.fragment
        or realm.username is not None
        or realm.password is not None
    ):
        raise ValueError("Registry Bearer realm is not the trusted GHCR token endpoint.")
    if attributes["service"] != GHCR_REGISTRY or attributes["scope"] != GHCR_PULL_SCOPE:
        raise ValueError("Registry Bearer challenge requested an unexpected service or scope.")
    return f"https://{GHCR_REGISTRY}/token?{urlencode({'service': GHCR_REGISTRY, 'scope': GHCR_PULL_SCOPE})}"


def _bearer_token(challenge: str, username: str, token: str, opener: Callable[..., object]) -> str:
    token_url = _trusted_token_url(challenge)
    result = _open(
        Request(token_url, headers={"Authorization": _basic_authorization(username, token)}), opener
    )
    if result.status != 200:
        raise ValueError(f"Registry token endpoint returned HTTP {result.status}.")
    try:
        payload = json.loads(result.body.decode("utf-8"))
        access_token = payload.get("token") or payload.get("access_token")
    except (UnicodeDecodeError, json.JSONDecodeError):
        access_token = None
    if not isinstance(access_token, str) or not access_token:
        raise ValueError("Registry token endpoint returned no usable bearer token.")
    return access_token


def check_registry_tag(
    registry: str,
    repository: str,
    tag: str,
    username: str,
    token: str,
    opener: Callable[..., object] = _open_without_redirects,
) -> RegistryResult:
    """Read one OCI manifest, permitting only explicit absence states."""

    if registry != GHCR_REGISTRY or repository != GHCR_REPOSITORY:
        return RegistryResult(RegistryDecision.ERROR, "unexpected_registry_or_repository")
    if not tag or not username or not token:
        return RegistryResult(RegistryDecision.ERROR, "missing_registry_credentials_or_identity")
    url = f"https://{registry}/v2/{repository}/manifests/{tag}"
    headers = {"Accept": OCI_MANIFEST_ACCEPT, "Authorization": _basic_authorization(username, token)}
    try:
        result = _open(Request(url, headers=headers), opener)
        if result.status == 401:
            challenge = result.headers.get("www-authenticate", "")
            bearer = _bearer_token(challenge, username, token, opener)
            result = _open(
                Request(url, headers={"Accept": OCI_MANIFEST_ACCEPT, "Authorization": f"Bearer {bearer}"}),
                opener,
            )
        return assess_manifest_response(result.status, result.body)
    except (TimeoutError, URLError, OSError, ValueError) as error:
        return RegistryResult(RegistryDecision.ERROR, f"registry_communication_error:{type(error).__name__}")


def _docker(*args: str, input_text: str | None = None) -> str:
    completed = subprocess.run(
        ["docker", *args], check=True, text=True, input=input_text, capture_output=True
    )
    return completed.stdout.strip()


def _label(image: str, name: str) -> str:
    return _docker("image", "inspect", "--format", f'{{{{ index .Config.Labels "{name}" }}}}', image)


def validate_image(
    image: str,
    expected_package_version: str,
    expected_container_version: str,
    expected_revision: str,
) -> None:
    """Validate the loaded release image which will be pushed without rebuilding."""

    configured_user = _docker("image", "inspect", "--format", "{{.Config.User}}", image)
    configured_uid = configured_user.split(":", 1)[0].strip().lower()
    if not configured_uid or configured_uid in {"root", "0"}:
        raise ValueError(f"Image is configured to run as root or has no configured user: {configured_user!r}")
    effective_uid = _docker("run", "--rm", image, "id", "-u")
    if not effective_uid.isdigit() or int(effective_uid) == 0:
        raise ValueError(f"Image runtime fell back to root (effective UID {effective_uid!r}).")

    expected_labels = {
        "org.opencontainers.image.source": "https://github.com/2002brian/nf-rna",
        "org.opencontainers.image.revision": expected_revision,
        "org.opencontainers.image.version": expected_container_version,
    }
    for name, expected in expected_labels.items():
        actual = _label(image, name)
        if actual.lower() in INVALID_REVISION_VALUES or actual != expected:
            raise ValueError(f"OCI label {name} is {actual!r}; expected {expected!r}.")

    cli_version = _docker("run", "--rm", image, "rnaseq", "--version")
    package_version = _docker(
        "run", "--rm", image, "python", "-c", "from rnaseq import __version__; print(__version__)"
    )
    if cli_version != expected_package_version or package_version != expected_package_version:
        raise ValueError(
            "Image package version mismatch: "
            f"CLI={cli_version!r}, package={package_version!r}, expected={expected_package_version!r}."
        )
    _docker("run", "--rm", image, "python", "-m", "rnaseq.workflow_support", "report", "--help")
    _docker(
        "run",
        "--rm",
        image,
        "python",
        "-c",
        "from rnaseq.workflow_assets import required_workflow_assets; "
        "missing=[name for name, path in required_workflow_assets().items() if not path.is_file()]; "
        "assert not missing, missing",
    )
    for asset in REQUIRED_R_ASSETS:
        _docker("run", "--rm", image, "test", "-f", f"/opt/nf-rna/r/{asset}")
    _docker(
        "run",
        "--rm",
        image,
        "sh",
        "-c",
        "command -v python && command -v Rscript && command -v ps",
    )
    _docker("run", "--rm", image, "python", "-c", "import rnaseq.models, rnaseq.workflow_support")
    _docker(
        "run",
        "--rm",
        image,
        "Rscript",
        "-e",
        "packages <- c(\"DESeq2\", \"tximport\", \"ggplot2\", \"pheatmap\", \"yaml\", "
        "\"jsonlite\", \"clusterProfiler\", \"AnnotationDbi\", \"org.Hs.eg.db\", \"org.Mm.eg.db\"); "
        "quit(status=if (all(vapply(packages, requireNamespace, logical(1), quietly=TRUE))) 0 else 1)",
    )


def _write_github_output(path: Path, values: Mapping[str, str]) -> None:
    with path.open("a", encoding="utf-8") as output:
        for name, value in values.items():
            output.write(f"{name}={value}\n")


def _release_info_command(args: argparse.Namespace) -> int:
    package_version = read_package_version(Path(args.package_version_file)) if args.package_version_file else None
    if args.resolve_tag and args.tag_revision:
        raise ValueError("Use either --resolve-tag or --tag-revision, not both.")
    tag_revision = resolve_tag_commit(args.tag) if args.resolve_tag else args.tag_revision
    identity = verify_release_identity(args.tag, package_version, tag_revision, args.checkout_revision)
    values = {
        "tag": identity.git_tag,
        "package_version": identity.package_version,
        "container_version": identity.container_version,
        "stable": str(identity.stable).lower(),
    }
    if tag_revision:
        values["sha"] = tag_revision
    if args.repository:
        values["image"] = f"ghcr.io/{args.repository}:{identity.container_version}"
    if args.github_output:
        _write_github_output(Path(args.github_output), values)
    else:
        print(json.dumps(values, sort_keys=True))
    return 0


def _registry_check_command(args: argparse.Namespace) -> int:
    token = os.environ.get(args.token_env, "")
    result = check_registry_tag(args.registry, args.repository, args.tag, args.username, token)
    print(json.dumps(asdict(result), sort_keys=True, default=str))
    return 0 if result.decision is RegistryDecision.ABSENT else 1


def _validate_image_command(args: argparse.Namespace) -> int:
    validate_image(
        args.image,
        args.expected_package_version,
        args.expected_container_version,
        args.expected_revision,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    release = commands.add_parser("release-info")
    release.add_argument("--tag", required=True)
    release.add_argument("--package-version-file")
    release.add_argument("--tag-revision")
    release.add_argument("--resolve-tag", action="store_true")
    release.add_argument("--checkout-revision")
    release.add_argument("--repository")
    release.add_argument("--github-output")
    release.set_defaults(func=_release_info_command)
    registry = commands.add_parser("registry-check")
    registry.add_argument("--registry", default="ghcr.io")
    registry.add_argument("--repository", required=True)
    registry.add_argument("--tag", required=True)
    registry.add_argument("--username", required=True)
    registry.add_argument("--token-env", default="GITHUB_TOKEN")
    registry.set_defaults(func=_registry_check_command)
    image = commands.add_parser("validate-image")
    image.add_argument("--image", required=True)
    image.add_argument("--expected-package-version", required=True)
    image.add_argument("--expected-container-version", required=True)
    image.add_argument("--expected-revision", required=True)
    image.set_defaults(func=_validate_image_command)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, subprocess.CalledProcessError) as error:
        print(f"GHCR release safeguard failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
