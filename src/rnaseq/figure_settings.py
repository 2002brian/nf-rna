"""Single resolver for figure-style settings.

Every producer of a figure configuration (direct L1/L2/GO/GSEA/KEGG runs, the
Nextflow config builders and the frozen downstream contract) goes through this
module, so the "absent means legacy" rule exists in exactly one place.

Resolution rules:

* ``project.yaml`` without ``figures``  -> ``legacy`` / ``general``.
* A frozen downstream contract without ``figures`` -> ``legacy`` / ``general``.
* An explicit selection is returned unchanged except that the alias ``bioinfo`` is
  normalised to the canonical preset ``standard`` (in every producer and every reader).
* One run has one setting: the frozen contract and the staged project snapshot
  must agree, otherwise :func:`resolve_run_figure_settings` raises.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Literal

from pydantic import ValidationError

from pydantic import field_validator

from rnaseq.models import PRESET_ALIASES, FiguresConfig, StrictModel

FIGURES_CONTRACT_KEY = "figures"


class ResolvedFigureSettings(StrictModel):
    """Fully resolved, serialisable figure settings carried to R."""

    # ``bioinfo`` is the alias of ``standard``: frozen contracts written before the rename
    # still parse and resolve to the canonical name (the frozen file itself is never rewritten).
    preset: Literal["standard", "bioinfo", "legacy"]
    profile: Literal["general", "nature", "plos", "elsevier"]
    # "explicit": selected in project.yaml; "legacy_default": nothing selected.
    source: Literal["explicit", "legacy_default"]

    @field_validator("preset")
    @classmethod
    def _canonical_preset(cls, value: str) -> str:
        return PRESET_ALIASES.get(value, value)


def _legacy() -> dict[str, str]:
    return {"preset": "legacy", "profile": "general", "source": "legacy_default"}


def _explicit(figures: FiguresConfig) -> dict[str, str]:
    return {"preset": figures.preset, "profile": figures.profile, "source": "explicit"}


def resolve_figure_settings(figures: FiguresConfig | Mapping[str, Any] | None) -> dict[str, str]:
    """Resolve a ``figures`` section (model, raw mapping, or absent) to settings."""

    if figures is None:
        return _legacy()
    if isinstance(figures, FiguresConfig):
        return _explicit(figures)
    if not isinstance(figures, Mapping):
        raise ValueError(f"invalid figures configuration: expected a mapping, got {type(figures).__name__}.")
    try:
        return _explicit(FiguresConfig.model_validate(dict(figures)))
    except ValidationError as exc:
        raise ValueError(f"invalid figures configuration: {exc}") from exc


def figure_settings_from_project(project: Mapping[str, Any]) -> dict[str, str]:
    """Resolve from a raw ``project.yaml`` mapping (``figures`` may be absent)."""

    return resolve_figure_settings(project.get(FIGURES_CONTRACT_KEY))


def figure_settings_from_contract(contract: Mapping[str, Any]) -> dict[str, str]:
    """Resolve from a frozen downstream contract; a missing key is legacy."""

    if FIGURES_CONTRACT_KEY not in contract or contract[FIGURES_CONTRACT_KEY] is None:
        return _legacy()
    try:
        return ResolvedFigureSettings.model_validate(contract[FIGURES_CONTRACT_KEY]).model_dump()
    except ValidationError as exc:
        raise ValueError(f"invalid figures in frozen downstream contract: {exc}") from exc


def resolve_run_figure_settings(contract: Mapping[str, Any], project: Mapping[str, Any]) -> dict[str, str]:
    """Return the run's single figure setting, rejecting a mixed-setting bundle.

    The frozen contract is authoritative (it is what the run was frozen with);
    the staged project snapshot must resolve to the same preset and profile.
    """

    from_contract = figure_settings_from_contract(contract)
    from_project = figure_settings_from_project(project)
    if (from_contract["preset"], from_contract["profile"]) != (from_project["preset"], from_project["profile"]):
        raise ValueError(
            "mixed figure settings in one run: the frozen downstream contract resolves to "
            f"{from_contract['preset']}/{from_contract['profile']} but the frozen project configuration resolves to "
            f"{from_project['preset']}/{from_project['profile']}."
        )
    return from_contract
