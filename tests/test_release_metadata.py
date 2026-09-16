"""Keep package, tag, and GHCR prerelease identities aligned."""

from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_ghcr_publication_accepts_pep440_prereleases_without_moving_latest():
    workflow = (ROOT / ".github" / "workflows" / "publish-ghcr.yml").read_text(encoding="utf-8")

    assert "^v[0-9]+\\.[0-9]+\\.[0-9]+((a|b|rc)[0-9]+)?$" in workflow
    assert 'if [[ "$version" =~ ^[0-9]+\\.[0-9]+\\.[0-9]+$ ]]; then' in workflow
    assert "type=raw,value=latest,enable=${{ steps.release.outputs.stable }}" in workflow
