"""L1 figures group samples by the same primary factor as Python (audit D-08).

Static check: the cloud container has no R.  The R-backed figure tests in the
main suite exercise the plotted output where Rscript is available.
"""

from __future__ import annotations

from pathlib import Path

from rnaseq.validators import condition_factor

L1 = Path(__file__).resolve().parents[2] / "src" / "rnaseq" / "r" / "l1_analysis.R"


def _r_choice(formula_terms: tuple[str, ...]) -> str:
    """Python model of the R expression that picks the L1 colour variable."""

    return "condition" if "condition" in formula_terms else formula_terms[-1]


def test_l1_colour_variable_prefers_condition_like_python():
    text = L1.read_text(encoding="utf-8")
    assert 'color_name <- if ("condition" %in% formula_terms) "condition" else tail(formula_terms, 1)' in text
    assert "color_name <- tail(all.vars(formula), 1)" not in text
    categorical = {"condition": "categorical", "batch": "categorical", "group": "categorical", "pair_id": "categorical"}
    for terms in [("condition",), ("batch", "condition"), ("condition", "batch"), ("pair_id", "condition"), ("batch", "group")]:
        assert _r_choice(terms) == condition_factor(terms, categorical)
