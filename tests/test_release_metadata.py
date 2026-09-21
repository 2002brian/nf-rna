"""Keep the release-only GHCR workflow simple and correctly ordered."""

from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_ghcr_publication_uses_standard_actions_for_a_stable_release():
    workflow = (ROOT / ".github" / "workflows" / "publish-ghcr.yml").read_text(encoding="utf-8")

    checkout = workflow.index("uses: actions/checkout@v6")
    verify = workflow.index("Verify release version and checkout")
    buildx = workflow.index("uses: docker/setup-buildx-action@v4")
    login = workflow.index("uses: docker/login-action@v4")
    build = workflow.index("uses: docker/build-push-action@v7")
    smoke = workflow.index("Smoke test published image")

    assert checkout < verify < buildx < login < build < smoke
    assert "workflow_dispatch:" not in workflow
    assert "^v[0-9]+\\.[0-9]+\\.[0-9]+$" in workflow
    assert "test \"$package_version\" = \"$expected_version\"" in workflow
    assert "test \"$checkout_revision\" = \"$release_revision\"" in workflow
    assert "contents: read" in workflow
    assert "packages: write" in workflow
    assert "platforms: linux/amd64" in workflow
    assert "NF_RNA_SOURCE_REVISION=${{ steps.release.outputs.sha }}" in workflow
    assert "ghcr.io/${{ github.repository }}:${{ steps.release.outputs.version }}" in workflow
    assert "ghcr.io/${{ github.repository }}:latest" in workflow
    assert "docker run --rm \"$IMAGE\" rnaseq --version" in workflow


def test_ghcr_publication_has_no_custom_registry_preflight_or_helper():
    workflow = (ROOT / ".github" / "workflows" / "publish-ghcr.yml").read_text(encoding="utf-8")

    assert "registry-check" not in workflow
    assert "tools/ghcr_release.py" not in workflow
    assert "docker manifest inspect" not in workflow
    assert "metadata-action" not in workflow
