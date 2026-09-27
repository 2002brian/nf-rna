"""Dependency invariants of the Conda environment that installs nf-rna."""

from __future__ import annotations

import ast
import re
import tomllib
from pathlib import Path

import yaml


ROOT = Path(__file__).parents[1]


def _normalized_distribution_name(specification: str) -> str:
    return re.split(r"[<>=!~;\[]", specification, maxsplit=1)[0].strip().lower().replace("_", "-")


def test_shared_execution_environment_covers_declared_python_runtime_dependencies():
    """environment.yml must supply every declared runtime dependency of nf-rna."""

    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    environment = yaml.safe_load((ROOT / "environment.yml").read_text(encoding="utf-8"))
    declared = {_normalized_distribution_name(item) for item in pyproject["project"]["dependencies"]}
    supplied = {_normalized_distribution_name(item) for item in environment["dependencies"] if isinstance(item, str)}

    assert declared <= supplied


def test_pydantic_execution_contract_requires_v2_apis():
    """Avoid a satisfiable-but-incompatible Pydantic v1 environment."""

    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    pydantic_requirement = next(
        item
        for item in pyproject["project"]["dependencies"]
        if _normalized_distribution_name(item) == "pydantic"
    )
    tree = ast.parse((ROOT / "src" / "rnaseq" / "models.py").read_text(encoding="utf-8"))
    imports = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module == "pydantic"
        for alias in node.names
    }

    assert {"ConfigDict", "field_validator", "model_validator"} <= imports
    assert ">=2.8" in pydantic_requirement and "<3" in pydantic_requirement
