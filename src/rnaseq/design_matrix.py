"""Dependency-free model-matrix dimensions shared by validation and the L2 launch guards."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence

_VARIABLE = re.compile(r"[A-Za-z_][A-Za-z0-9_.]*")


def formula_terms(formula: str) -> tuple[str, ...]:
    """Return the additive variables of a restricted formula such as ``~ batch + condition``."""

    return tuple(_VARIABLE.findall(formula.removeprefix("~").strip()))


def model_coefficient_count(
    rows: Sequence[Mapping[str, str]], terms: Iterable[str], variable_types: Mapping[str, str]
) -> int | None:
    """Count the columns of the expanded additive model matrix: intercept, one slope per
    continuous variable, and one column per non-baseline level of each categorical variable.

    Returns ``None`` when a term is untyped or absent from the rows, so callers never
    block on a design they cannot interpret (for example an older project that
    declared no variable types).
    """

    coefficients = 1
    for term in terms:
        kind = variable_types.get(term)
        if kind is None or any(term not in row for row in rows):
            return None
        if kind == "continuous":
            coefficients += 1
        else:
            coefficients += max(len({row[term] for row in rows}) - 1, 0)
    return coefficients


def residual_degrees_of_freedom_message(
    rows: Sequence[Mapping[str, str]], formula: str, variable_types: Mapping[str, str]
) -> str | None:
    """Return a user-facing message when the fitted design leaves no residual degrees of freedom."""

    coefficients = model_coefficient_count(rows, formula_terms(formula), variable_types)
    if coefficients is None or not rows or len(rows) > coefficients:
        return None
    return (
        f"The design {formula} has no residual degrees of freedom: {len(rows)} sample(s) for {coefficients} "
        "model coefficient(s) after expanding categorical variables, so DESeq2 cannot estimate dispersion."
    )
