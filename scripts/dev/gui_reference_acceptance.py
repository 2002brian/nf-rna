"""New Project reference selection in headless Chrome (v1.7.0 Human HISAT2 reference fix).

Usage (the GUI must already be running; it is not started here):

    python scripts/dev/gui_reference_acceptance.py --base URL --parent DIR --shots DIR --evidence FILE

Uses this workstation's real reference registry (read only). Checks Human/Mouse x Salmon/HISAT2,
the species and backend switches, then creates one Human + HISAT2 FASTQ project scaffold in
--parent (an allowed root of the running GUI) and opens Validate & Run for it. Nothing is run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml
from playwright.sync_api import Page, sync_playwright

sys.path.insert(0, str(Path(__file__).parent))
import gui_browser_acceptance as gui  # noqa: E402
import gui_workspace_acceptance as ws  # noqa: E402

HUMAN = "/home/brian/reference/homo_sapiens/ensembl_116_GRCh38"
MOUSE = "/home/brian/reference/mus_musculus/ensembl_116_GRCm39"


def reference_value(page: Page) -> str:
    return page.get_by_label("Reference", exact=True).input_value()  # the selectbox shows its choice in the combobox input


def state(page: Page, species: str, method: str) -> str:
    ws.pick(page, "Species", species)
    ws.pick(page, "Quantification", method)
    return reference_value(page)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--shots", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    gui.BASE = args.base
    args.shots.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(executable_path=gui.CHROME, headless=True)
        page = browser.new_context(viewport={"width": 1500, "height": 1100}).new_page()
        gui.goto(page, "new", project=None)
        page.get_by_text("fastq — Reads requiring upstream processing", exact=True).click()
        gui.settle(page)
        for species, method, root, other in (
            ("Homo sapiens", "Salmon (pseudoalignment)", HUMAN, MOUSE), ("Homo sapiens", "HISAT2 + featureCounts (genome alignment)", HUMAN, MOUSE),
            ("Mus musculus", "Salmon (pseudoalignment)", MOUSE, HUMAN), ("Mus musculus", "HISAT2 + featureCounts (genome alignment)", MOUSE, HUMAN),
        ):
            shown = state(page, species, method)
            ws.expect(f"reference: {species} + {method} shows its own registered reference", root in shown and other not in shown, shown)
        shown = state(page, "Homo sapiens", "HISAT2 + featureCounts (genome alignment)")
        ws.expect("reference: Mouse -> Human switch shows the Human reference", HUMAN in shown, shown)
        shown = state(page, "Mus musculus", "HISAT2 + featureCounts (genome alignment)")
        ws.expect("reference: Human -> Mouse switch shows the Mouse reference", MOUSE in shown, shown)
        ws.pick(page, "Species", "Homo sapiens")
        ws.pick(page, "Quantification", "Salmon (pseudoalignment)")
        ws.pick(page, "Reference", "iGenomes")
        shown = state(page, "Homo sapiens", "HISAT2 + featureCounts (genome alignment)")
        ws.expect("reference: Salmon -> HISAT2 does not keep iGenomes", "iGenomes" not in shown and HUMAN in shown, shown)
        ws.expect("reference: no 'not prepared' notice for a registered Human HISAT2 reference", "not prepared" not in ws.main_text(page), "")
        gui.shot(page, args.shots, "R1_human_hisat2_reference")

        name = "human-hisat2-ref-check"
        gui.fill(page, "Project name", name)
        gui.fill(page, "Location (parent folder)", str(args.parent))
        ws.pick(page, "Strandedness", "reverse")
        gui.click(page, "Create project")
        created = args.parent / name
        config = yaml.safe_load((created / "project.yaml").read_text()) if (created / "project.yaml").is_file() else {}
        ws.expect("reference: Human + HISAT2 project created with the Human root",
                  config.get("reference", {}).get("root") == HUMAN and config.get("organism", {}).get("species") == "Homo sapiens", config.get("reference"))
        gui.goto(page, "run", project=str(created))
        text = ws.main_text(page)
        reference_errors = [line for line in text.splitlines() if "reference" in line.lower() and ("must" in line.lower() or "requires" in line.lower() or "not match" in line.lower())]
        from rnaseq.validators import validate_project

        report = validate_project(created)
        loaded = report.local_reference
        ws.expect("reference: the project validator loads the Human HISAT2 reference (runtime-ready) with no reference error",
                  loaded is not None and loaded.species == "Homo sapiens" and loaded.hisat2_runtime_ready
                  and not [issue.code for issue in report.errors if "reference" in issue.code or "hisat2" in issue.code],
                  {"assembly": loaded.assembly_identity if loaded else None, "errors": [issue.code for issue in report.errors]})
        ws.expect("reference: no reference error among the findings", not reference_errors, reference_errors)
        gui.shot(page, args.shots, "R2_validate_human_hisat2_project")
        browser.close()
    args.evidence.write_text(json.dumps(ws.EVIDENCE, indent=2) + "\n")
    failed = [key for key, item in ws.EVIDENCE.items() if not item["pass"]]
    print(f"{len(ws.EVIDENCE) - len(failed)}/{len(ws.EVIDENCE)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
