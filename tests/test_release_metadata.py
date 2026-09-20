"""Keep the workflow wired to the behavioral GHCR safeguard helper."""

from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_ghcr_publication_uses_canonical_release_helper_and_stable_only_latest():
    workflow = (ROOT / ".github" / "workflows" / "publish-ghcr.yml").read_text(encoding="utf-8")

    checkout = workflow.index("uses: actions/checkout@v6")
    helper = workflow.index("python3 tools/ghcr_release.py")
    build = workflow.index("uses: docker/build-push-action@v7")

    assert workflow.count("python3 tools/ghcr_release.py release-info") == 1
    assert checkout < helper < build
    assert "tools/ghcr_release.py" not in workflow[:checkout]
    assert "context: ." in workflow[build:]
    assert "--resolve-tag" in workflow
    assert "python3 tools/ghcr_release.py validate-image" in workflow
    assert "python3 tools/ghcr_release.py registry-check" in workflow
    assert "container_version" in workflow
    assert "type=raw,value=latest,enable=${{ steps.release.outputs.stable }}" in workflow


def test_ghcr_publication_is_release_only_and_checks_registry_before_login_or_push():
    workflow = (ROOT / ".github" / "workflows" / "publish-ghcr.yml").read_text(encoding="utf-8")

    assert "workflow_dispatch:" not in workflow
    assert "docker manifest inspect" not in workflow
    assert workflow.index("Refuse to replace an existing release image") < workflow.index("Set up Docker Buildx")
    assert workflow.index("Refuse to replace an existing release image") < workflow.index("Log in to GitHub Container Registry")
    assert workflow.index("Log in to GitHub Container Registry") < workflow.index("Publish validated image tags")
    publish = workflow[workflow.index("- name: Publish validated image tags") :]
    assert publish.index('docker push "$IMAGE"') < publish.index('docker push "$LATEST_IMAGE"')
    assert 'if [[ "$STABLE" == "true" ]]; then' in publish
