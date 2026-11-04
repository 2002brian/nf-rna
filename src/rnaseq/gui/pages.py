"""Streamlit pages: Projects, Datasets, Project, New Project, Input & Metadata, Contrasts, Run, Results, Reports & Downloads.

Pages read every state from the project and run files on each rerun;
``st.session_state`` holds only unsaved UI input (the metadata being edited,
the contrasts being composed, the selected project).
"""

from __future__ import annotations

import shlex
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

from rnaseq.errors import ExecutionPreflightError, ProjectConfigError, ProjectCreationError
from rnaseq.gui import datasets, figure_export, jobs, results, services, workspace
from rnaseq.gui.config import GuiConfig, PathNotAllowed, ensure_data_layout, gui_state_dir, load_config
from rnaseq.references import enable_checksum_memo
from rnaseq.models import DEFAULT_L2_ENRICHMENT, PIPELINE_VERSION, DesignType, FastqPreprocessing, InputType, Preset, SequencingLayout, describe_enrichment
from rnaseq.run_status import _Cache, build_status, discover_runs, render_status
from rnaseq import wizard

PROJECT_KEY = "project"

STRANDEDNESS_HELP = """\
**HISAT2 + featureCounts needs the library strandedness; it is never guessed.**

- **reverse** — dUTP / "directional" stranded kits, read 1 antisense: Illumina TruSeq Stranded mRNA and
  Total RNA, Illumina Stranded mRNA Prep, NEBNext Ultra II Directional. Most stranded Illumina libraries are this.
- **forward** — read 1 sense: for example Lexogen QuantSeq 3′ mRNA FWD, and ligation-based kits documented as sense-stranded.
- **unstranded** — non-stranded kits, for example Illumina TruSeq RNA v2 (non-stranded), NEBNext Ultra II (non-directional), SMART-Seq v4.

Not sure? Check the kit name in the sequencing provider's report or the kit manual, ask the sequencing core,
or run the project once with Salmon (strandedness `auto` infers it), or RSeQC `infer_experiment.py` on a few BAM files.
"""


BACKEND_HELP = {
    "salmon": "Salmon: transcript-level quantification by pseudoalignment (nf-core/rnaseq), summarised to genes with tximport. "
              "Strandedness `auto` lets Salmon infer it from the reads.",
    "hisat2_featurecounts": "HISAT2 + featureCounts: splice-aware genome alignment, then gene-level read counting on the annotation. "
                            "featureCounts needs the strandedness stated explicitly.",
}

DESIGN_HELP = {
    "two_group": "Unpaired: independent samples in two conditions (for example different animals per group). The model is `~ condition`.",
    "paired_two_group": "Paired: each pairing unit (patient, donor, cell line …) has one sample of each of the two conditions, "
                        "for example tumour and adjacent normal of the same patient. The model compares the conditions within each pair: `~ <pairing column> + condition`.",
    "multi_group": "Three or more condition levels with independent samples; contrasts between levels are chosen on the Contrasts page.",
}


# --------------------------------------------------------------------- shared context

def gui_config() -> GuiConfig:
    # Every rerun re-validates the project, which re-checks the reference assets'
    # SHA-256 (several GB).  Within this long-lived process, reuse the digest of
    # an unchanged large file (identity = device, inode, size, mtime, ctime);
    # `rnaseq run` and the CLI still hash in full.  See references.sha256_file.
    enable_checksum_memo()
    config = load_config()
    if not st.session_state.get("_layout_checked"):
        st.session_state["_layout_messages"] = [*config.notes, *ensure_data_layout(config)]
        st.session_state["_layout_checked"] = True
    return config


def current_project() -> Path | None:
    """The selected project, from the session or (on a reload or bookmark) the URL.

    A project that arrives by URL is copied into the session, and the URL is
    restored after a page switch (``st.switch_page`` and page links drop the
    query string), so changing pages never loses the project.
    """

    value = st.session_state.get(PROJECT_KEY) or st.query_params.get(PROJECT_KEY)
    if not value:
        return None
    path = allowed_project(load_config(), value)
    if path is None:
        # A URL (or a stale session value) naming a folder outside the allowed
        # roots is dropped; it is never opened, read or written.
        st.session_state.pop(PROJECT_KEY, None)
        if PROJECT_KEY in st.query_params:
            st.query_params.pop(PROJECT_KEY, None)
        st.session_state["_project_refused"] = f"The project `{value}` was not opened: it is not a project inside an available data folder (allowed roots in the GUI settings)."
        return None
    st.session_state[PROJECT_KEY] = str(path)
    if st.query_params.get(PROJECT_KEY) != str(path):
        st.query_params[PROJECT_KEY] = str(path)
    return path


def allowed_project(config: GuiConfig, value: str | Path) -> Path | None:
    """``value`` resolved, when it is a project folder inside an available allowed root (else None)."""

    try:
        path = config.resolve_allowed(value)
    except PathNotAllowed:
        return None
    return path if (path / "project.yaml").is_file() else None


def select_project(path: Path | None) -> None:
    if path is None:
        st.session_state.pop(PROJECT_KEY, None)
        st.query_params.pop(PROJECT_KEY, None)
        return
    st.session_state[PROJECT_KEY] = str(path)
    st.query_params[PROJECT_KEY] = str(path)


def require_project() -> Path | None:
    project = current_project()
    if project is None:
        refused = st.session_state.pop("_project_refused", None)
        if refused:
            st.error(refused)
        st.info("Select a project in the sidebar, or create one on the New Project page.")
    return project


VALIDATING = ("Validating the project … The first check of a reference in this GUI session verifies its SHA-256 "
              "checksums (several GB) and can take 5–10 s; later checks reuse them while the files are unchanged.")


def validated(project: Path):
    """``validate_project`` (read-only); a note appears if it takes longer than half a second."""

    from rnaseq.validators import validate_project

    with st.spinner(VALIDATING):
        return validate_project(project)


def workflow_bar(project: Path, report, here: str) -> None:
    """Where this page sits in the workflow, what is done, and what to do next (read from the project files)."""

    if report.config is None or report.loaded is None:
        return
    steps = workspace.workflow_steps(report, workspace.project_runs(project, _Cache()))
    parts = []
    for step in steps:
        name = f"**{step.name}**" if STEP_PAGES.get(step.name) == here else step.name
        parts.append(f"{STEP_ICONS.get(step.status, '')} {name}")
    box = st.container(border=True)
    box.markdown(" → ".join(parts), help="✅ done · ⬜ to do · ⚠️ needs attention · ➖ not needed for this analysis scope. "
                 "Bold steps are done on this page. The state is read from the project files on every visit.")
    upcoming = workspace.next_step(steps)
    if upcoming is None:
        box.caption("Every step is complete.")
        return
    target = STEP_PAGES.get(upcoming.name)
    text = f"Next: **{upcoming.name}** — {upcoming.detail}"
    if target == here:
        box.caption(text + " (on this page).")
    elif target in PAGES:
        cols = box.columns([3, 1])
        cols[0].caption(text + ".")
        cols[1].page_link(PAGES[target], label=f"Go to {PAGES[target].title} →")
    else:
        box.caption(text + ".")


def sidebar(config: GuiConfig) -> None:
    projects = services.discover_projects(config)
    current = current_project()
    options = [str(path) for path in projects]
    if current is not None and str(current) not in options:
        options.insert(0, str(current))
    with st.sidebar:
        st.caption(f"nf-rna {PIPELINE_VERSION}")
        refused = st.session_state.pop("_project_refused", None)
        if refused:
            st.error(refused)
        if options:
            # The widget follows the current project, which other pages can change.
            st.session_state["sidebar_project"] = str(current) if current is not None else None
            st.selectbox("Project", options, index=None, placeholder="Choose a project", key="sidebar_project",
                         format_func=lambda value: Path(value).name,
                         on_change=lambda: select_project(Path(st.session_state["sidebar_project"])) if st.session_state.get("sidebar_project") else None)
        with st.expander("Open project by path"):
            typed = st.text_input("Project directory", key="sidebar_open_path")
            if st.button("Open", key="sidebar_open"):
                try:
                    path = config.resolve_allowed(typed)
                    if not (path / "project.yaml").is_file():
                        raise PathNotAllowed(f"{path} has no project.yaml.")
                    select_project(path)
                    st.rerun()
                except PathNotAllowed as exc:
                    st.error(str(exc))
        if current is not None:
            st.caption(f"Current: `{current}`")
            summary = _sidebar_summary(str(current))
            if summary:
                st.caption(summary)


def _sidebar_summary(project: str) -> str:
    summary = workspace.project_summary(Path(project))
    if summary.error:
        return "project.yaml cannot be read"
    latest = summary.latest
    state = f"latest run {latest.execution}, QC {latest.qc.headline()}" if latest else "no run yet"
    return f"{summary.scope} · {summary.input_type} · {summary.samples} samples · {state}"


def remembered_choice(label: str, options: list[str], key: str, **kwargs) -> str:
    """A selectbox whose choice survives visits to other projects and pages.

    Streamlit drops the state of a widget that is not drawn in a rerun, so the
    choice is also kept under ``_keep::<key>`` and restored when the widget is
    drawn again (if it is still one of the options).
    """

    keep = f"_keep::{key}"
    if st.session_state.get(key) not in options:
        st.session_state.pop(key, None)
        if st.session_state.get(keep) in options:
            st.session_state[key] = st.session_state[keep]
    value = st.selectbox(label, options, index=0, key=key, **kwargs)
    st.session_state[keep] = value
    return value


def remembered_text(label: str, key: str, **kwargs) -> str:
    """A text input whose value survives visits to other projects and pages (see :func:`remembered_choice`)."""

    keep = f"_keep::{key}"
    if key not in st.session_state and keep in st.session_state:
        st.session_state[key] = st.session_state[keep]
    value = st.text_input(label, key=key, **kwargs)
    st.session_state[keep] = value
    return value


def remember(key: str, value: str) -> None:
    """Preselect ``value`` for the :func:`remembered_choice` widget ``key``."""

    st.session_state[key] = value
    st.session_state[f"_keep::{key}"] = value


def go_to(page: str) -> None:
    """Switch page when running inside the navigation (not in a single-page AppTest)."""

    if page in PAGES:
        st.switch_page(PAGES[page])


def path_picker(config: GuiConfig, label: str, key: str, *, want: str, suffixes: tuple[str, ...] = (), start: Path | None = None) -> Path | None:
    """A text field plus a directory listing; only paths inside the allowed roots are returned."""

    roots = config.available_roots()
    if not roots:
        st.error("No allowed root is available; check the GUI settings and mounted disks.")
        return None
    text_key, cwd_key, error_key = f"{key}_text", f"{key}_cwd", f"{key}_pick_error"
    resolved_roots = [root.resolve() for root in roots]
    if cwd_key not in st.session_state:
        st.session_state[cwd_key] = str(start if start is not None and start.is_dir() else roots[0])
    # The browsed folder is always a resolved path inside an available root; a
    # stale or forged value (an unmounted disk, a link) falls back to the first root.
    try:
        cwd = config.resolve_allowed(st.session_state[cwd_key])
        if not cwd.is_dir():
            raise PathNotAllowed(f"{cwd} is not a folder.")
    except PathNotAllowed:
        cwd = resolved_roots[0]
    st.session_state[cwd_key] = str(cwd)

    # Buttons change other widgets' values only in callbacks, which run before
    # the widgets are created in the next rerun (Streamlit forbids changing the
    # value of a widget that was already drawn in the current run).
    def browse_to(target: Path) -> None:
        try:
            resolved = config.resolve_allowed(target)
        except PathNotAllowed as exc:
            st.session_state[error_key] = f"Not opened: {exc}"
            return
        st.session_state[cwd_key] = str(resolved)
        st.session_state.pop(error_key, None)

    def use(target: Path, kind: str) -> None:
        try:
            resolved = config.resolve_allowed(target)
            if (kind == "dir" and not resolved.is_dir()) or (kind == "file" and not resolved.is_file()):
                raise PathNotAllowed(f"{resolved} is not a {'folder' if kind == 'dir' else 'file'}.")
        except PathNotAllowed as exc:
            st.session_state[error_key] = f"Not selected: {exc}"
            return
        st.session_state[text_key] = str(resolved)
        st.session_state.pop(error_key, None)

    def use_file() -> None:
        name = st.session_state.get(f"{key}_file")
        if name:
            use(cwd / name, "file")

    def open_sub() -> None:
        name = st.session_state.get(f"{key}_sub")
        if name:
            browse_to(cwd / name)

    typed = st.text_input(label, key=text_key, placeholder=f"A {'directory' if want == 'dir' else 'file'} inside {', '.join(map(str, roots))}")
    if st.session_state.get(error_key):
        st.error(st.session_state[error_key])
    with st.expander("Browse", expanded=not typed):
        top = st.columns([1, 3])
        root_choice = top[0].selectbox("Root", [str(root) for root in roots], key=f"{key}_root")
        top[0].button("Go to root", key=f"{key}_goroot", on_click=browse_to, args=(Path(root_choice),))
        top[1].markdown(f"Folder: `{cwd}`")
        try:
            entries = sorted(cwd.iterdir(), key=lambda item: item.name.lower())
        except OSError as exc:
            st.error(f"Cannot list {cwd}: {exc}")
            entries = []
        visible = [entry for entry in entries if not entry.name.startswith(".")]
        inside = []
        for entry in visible:
            try:
                config.resolve_allowed(entry)
            except PathNotAllowed:
                continue
            inside.append(entry)
        directories = [entry for entry in inside if entry.is_dir()]
        files = [entry for entry in inside if entry.is_file() and (not suffixes or entry.name.endswith(suffixes))]
        if len(inside) < len(visible):
            st.caption(f"{len(visible) - len(inside)} link(s) pointing outside the allowed data folders are not shown.")
        row = st.columns([3, 1, 1])
        sub = row[0].selectbox("Subfolders", [entry.name for entry in directories], index=None, placeholder="(none selected)", key=f"{key}_sub")
        row[1].button("Open folder", key=f"{key}_open", disabled=sub is None, on_click=open_sub)
        row[2].button("Up", key=f"{key}_up", disabled=cwd in resolved_roots, on_click=browse_to, args=(cwd.parent,))
        if want == "dir":
            st.button("Use this folder", key=f"{key}_usedir", on_click=use, args=(cwd, "dir"))
        else:
            row = st.columns([3, 1])
            chosen = row[0].selectbox("Files", [entry.name for entry in files], index=None, placeholder="(none selected)", key=f"{key}_file")
            row[1].button("Use this file", key=f"{key}_usefile", disabled=chosen is None, on_click=use_file)
    if not typed:
        return None
    try:
        path = config.resolve_allowed(typed)
    except PathNotAllowed as exc:
        st.error(str(exc))
        return None
    if want == "dir" and not path.is_dir():
        st.error(f"Not a directory: {path}")
        return None
    if want == "file" and not path.is_file():
        st.error(f"Not a file: {path}")
        return None
    return path


def show_issues(issues, *, title: str | None = None) -> None:
    for issue in issues:
        (st.error if issue.severity.value == "ERROR" else st.warning)(issue.message)


# --------------------------------------------------------------------- A. Dashboard

@st.cache_data(ttl=600, show_spinner="Running doctor checks ...")
def _doctor() -> list[tuple[str, str, str]]:
    from rnaseq.execution import doctor_checks

    return [(item.name, item.verdict, item.detail) for item in doctor_checks(None)]


def _project_table_row(item: workspace.ProjectSummary) -> dict[str, object]:
    latest = item.latest
    # Status columns first: the table is wider than most screens.
    return {
        "project": item.project_id, "scope": item.scope, "samples": item.samples, "contrasts": ", ".join(item.contrasts) or "—",
        "last run": latest.started[:16].replace("T", " ") if latest else "—",
        "execution": latest.execution if latest else "—", "scientific QC": latest.qc.headline() if latest else "—",
        "runs": len(item.runs), "species": item.species, "input": item.input_type, "quantification": item.quantification,
        "reference": item.reference, "path": str(item.path), "problem": item.error or "",
    }


def dashboard_page() -> None:
    config = gui_config()
    st.title("Research projects")
    st.caption(f"nf-rna {PIPELINE_VERSION} · settings: `{config.source}`")
    for message in st.session_state.get("_layout_messages", []):
        st.info(message)
    active = jobs.active_launch()
    if active is not None:
        st.warning(f"Active GUI run: case **{active.case_id}** of `{active.project}` (submitted {active.launched_at}). Open the Run page of that project to follow it.")

    for item in config.root_statuses():
        if not item.available:
            st.error(f"Data folder `{item.path}` is unavailable: {item.message} Its projects and datasets are not shown.")

    st.subheader("Research projects")
    cache = _Cache()
    projects = services.discover_projects(config)
    summaries = [workspace.project_summary(path, cache) for path in projects]
    if summaries:
        st.dataframe(pd.DataFrame([_project_table_row(item) for item in summaries]), hide_index=True, width="stretch",
                     column_config={"path": st.column_config.TextColumn("path", width="small")})
        cols = st.columns([3, 1])
        chosen = cols[0].selectbox("Open project", [str(item.path) for item in summaries], index=None, format_func=lambda value: Path(value).name,
                                   placeholder="Choose a project to open", key="dash_open")
        cols[1].write("")
        if cols[1].button("Open", key="dash_open_btn", disabled=chosen is None, type="primary"):
            select_project(Path(chosen))
            go_to("project")
            st.rerun()
    else:
        st.write("No projects found under the allowed roots. Create one on **New Project**, or open one by path in the sidebar.")

    st.subheader("Recent runs")
    rows = []
    for item in summaries:
        for run in item.runs[:15]:
            rows.append({"started": run.started, "project": item.project_id, "run": run.label, "level": run.level or "",
                         "execution": run.execution, "scientific QC": run.qc.headline(), "note": run.note or ""})
    rows.sort(key=lambda row: row["started"], reverse=True)
    if rows:
        st.dataframe(pd.DataFrame(rows[:15]), hide_index=True, width="stretch")
    else:
        st.write("No recorded runs.")

    st.subheader("Server")
    st.dataframe(pd.DataFrame([{"data folder (allowed root)": str(item.path), "status": "available" if item.available else "UNAVAILABLE", "detail": item.message} for item in config.root_statuses()]), hide_index=True, width="stretch")
    poppler, message = figure_export.poppler_status()
    (st.caption if poppler else st.info)(message)
    left, right = st.columns(2)
    with left:
        with st.expander("Runtime checks (rnaseq doctor)"):
            if st.button("Re-run checks", key="doctor_refresh"):
                _doctor.clear()
            checks = _doctor()
            st.dataframe(pd.DataFrame(checks, columns=["check", "verdict", "detail"]), hide_index=True, width="stretch")
    with right, st.expander("Registered references"):
        references, error = services.registered_references()
        if error:
            st.error(error)
        if references:
            st.dataframe(pd.DataFrame([{"species": ref.species, "provider": ref.provider, "release": ref.release, "assembly": ref.assembly_identity, "Salmon": "yes" if ref.salmon_index else "no", "HISAT2": "yes" if ref.hisat2_index else "no", "root": str(ref.root)} for ref in references]), hide_index=True, width="stretch")
        elif not error:
            st.write("No registered references.")


# --------------------------------------------------------------------- Project workspace

STEP_ICONS = {"done": "✅", "todo": "⬜", "problem": "⚠️", "not needed": "➖"}
STEP_PAGES = {"Dataset": "inputs", "Project": "new", "Metadata": "inputs", "Contrasts": "contrasts", "Validate": "run", "Run": "run", "Results": "results"}


def run_history_table(runs: list[workspace.RunRow]) -> pd.DataFrame:
    return pd.DataFrame([{
        "run": run.label, "level": run.level or "", "execution": f"{workspace.EXECUTION_ICONS.get(run.execution, '')} {run.execution}",
        "scientific QC": f"{workspace.QC_ICONS.get(run.qc.level, '')} {run.qc.headline()}",
        "started": run.started[:19].replace("T", " "), "completed": run.completed[:19].replace("T", " "),
        "results": "yes" if run.has_results else "", "report": "yes" if run.has_report else "", "note": run.note or "",
    } for run in runs])


def show_run_qc(qc: workspace.RunQc) -> None:
    """The recorded scientific QC verdicts of a run (display only; the execution status is separate)."""

    text = {
        "PASS": "Scientific QC: **PASS** — every recorded check passed.",
        "WARNING": f"Scientific QC: **WARNING** — {len(qc.warnings)} recorded check(s) need attention. The run's execution status is unchanged.",
        "INCONCLUSIVE": "Scientific QC: **INCONCLUSIVE** — at least one check had too little evidence for a verdict.",
        "NOT_EVALUATED": "Scientific QC: **NOT EVALUATED** — this run recorded no QC verdict (it did not reach the analysis, or it predates the check).",
    }[qc.level]
    {"PASS": st.success, "WARNING": st.warning, "INCONCLUSIVE": st.info, "NOT_EVALUATED": st.caption}[qc.level](text)
    if qc.items:
        st.dataframe(pd.DataFrame([{"check": item.check, "QC": f"{workspace.QC_ICONS.get(item.level, '')} {item.level}", "recorded": item.recorded, "detail": item.detail} for item in qc.items]),
                     hide_index=True, width="stretch", column_config={"detail": st.column_config.TextColumn("detail", width="large")})


def open_run(project: Path, run: workspace.RunRow, page: str) -> None:
    """Preselect ``run`` on the Results or Reports page and go there."""

    remember(f"{page}_run::{project}", run.label)
    go_to(page)


def project_page() -> None:
    gui_config()
    project = require_project()
    if project is None:
        return
    summary = workspace.project_summary(project)
    st.title(summary.project_id)
    st.caption(f"`{project}`")
    if summary.error:
        st.error(f"project.yaml cannot be read: {summary.error}")
        return
    cols = st.columns(4)
    cols[0].metric("Species", summary.species)
    cols[1].metric("Input", summary.input_type.replace("_", " "))
    cols[2].metric("Samples", summary.samples)
    cols[3].metric("Contrasts", len(summary.contrasts))
    facts_table({
        "quantification": summary.quantification, "reference": summary.reference,
        "layout / strandedness": f"{summary.layout} / {summary.strandedness}" if summary.input_type == InputType.FASTQ.value else None,
        "analysis": summary.analysis_settings(), "contrasts": ", ".join(summary.contrasts) or None,
    })

    report = validated(project)
    st.subheader("Workflow")
    if report.config is None:
        show_issues(report.errors)
        return
    steps = workspace.workflow_steps(report, list(summary.runs))
    cols = st.columns(len(steps))
    for col, step in zip(cols, steps):
        col.markdown(f"{STEP_ICONS.get(step.status, '')} **{step.name}**  \n{step.detail}")
    upcoming = workspace.next_step(steps)
    if upcoming is None:
        st.success("Every step is complete. Open the results of a run below.")
    else:
        st.info(f"Next: **{upcoming.name}** — {upcoming.detail}.")
        if upcoming.issues:
            show_issues(upcoming.issues[:8])
        target = STEP_PAGES.get(upcoming.name)
        if target in PAGES:
            st.page_link(PAGES[target], label=f"Go to {PAGES[target].title} →")

    with st.expander("Data locations"):
        assert report.loaded is not None
        origin = workspace.input_origin(project)
        facts_table({
            "project input": str(report.loaded.input_path),
            "original FASTQ location": ", ".join(origin) if origin else None,
            "metadata": str(report.loaded.metadata_path), "contrasts": str(report.loaded.contrasts_path),
            "runs": str(project / "runs"),
        })
        if report.config.input.type is InputType.RAW_COUNTS:
            st.caption("The count matrix was copied into the project; the importer does not record its original file.")

    st.subheader("Analysis history")
    runs = list(summary.runs)
    if not runs:
        st.write("No run yet.")
        if "run" in PAGES:
            st.page_link(PAGES["run"], label="Validate and run →")
        return
    st.dataframe(run_history_table(runs), hide_index=True, width="stretch")
    labels = [run.label for run in runs]
    chosen = remembered_choice("Run", labels, f"project_run::{project}")
    run = runs[labels.index(chosen)]
    show_run_qc(run.qc)
    cols = st.columns(3)
    if cols[0].button("Open results", key=f"project_open_results::{project}", disabled=not run.has_results):
        open_run(project, run, "results")
    if cols[1].button("Open report & downloads", key=f"project_open_reports::{project}", disabled=not run.has_results):
        open_run(project, run, "reports")
    if cols[2].button("Run details and logs", key=f"project_open_run::{project}"):
        remember(f"run_selected::{project}", run.label)
        go_to("run")


# --------------------------------------------------------------------- Datasets

@st.cache_data(ttl=60, max_entries=16, show_spinner="Listing the folder ...")
def _scan(folder: str, roots: tuple[str, ...]) -> datasets.DatasetScan:
    return datasets.scan_dataset(load_config(), folder)


@st.cache_data(max_entries=8, show_spinner="Checking the count matrix ...")
def _recognise_counts(path: str, mtime: float):
    return datasets.recognise_counts(Path(path))


def copy_plan_notice(plan: datasets.CopyPlan) -> None:
    free = datasets.human_size(plan.free_bytes) if plan.free_bytes is not None else "unknown"
    text = (f"Importing copies {plan.files} FASTQ file(s), {datasets.human_size(plan.bytes)}, into `{plan.target}` "
            f"(free space there: {free}). The original files stay where they are and are not changed.")
    (st.info if plan.fits else st.error)(text if plan.fits else text + " There is not enough free space.")


def import_fastq(config: GuiConfig, project: Path, records, sheet_text: str | None = None, sheet_path: Path | None = None) -> None:
    """Copy FASTQ files into the project with the existing importer (``rnaseq new --fastq-samplesheet``).

    Refused when a file resolves outside the allowed roots or the target has no room.
    """

    from rnaseq.project import import_fastq_samplesheet, load_project

    outside = datasets.outside_allowed(config, datasets.fastq_paths(records))
    if outside:
        st.error("These FASTQ files are outside the allowed data folders (for example through a symbolic link) and are not imported: " + ", ".join(outside[:5]))
        return
    try:
        config.resolve_allowed(project)
    except PathNotAllowed as exc:
        st.error(f"The project is not inside an available data folder; nothing is copied. {exc}")
        return
    plan = datasets.fastq_copy_plan(records, load_project(project).input_path)
    if not plan.fits:
        st.error(f"Not enough free space in `{plan.target}` for {datasets.human_size(plan.bytes)}; nothing is copied.")
        return
    if sheet_path is not None:
        _import_samplesheet(project, sheet_path, import_fastq_samplesheet, plan)
        return
    assert sheet_text is not None
    scratch = gui_state_dir() / "samplesheets" / f"{project.name}-{datetime.now():%Y%m%d-%H%M%S}.csv"
    scratch.parent.mkdir(parents=True, exist_ok=True)
    scratch.write_text(sheet_text, encoding="utf-8")
    try:
        _import_samplesheet(project, scratch, import_fastq_samplesheet, plan)
    finally:
        scratch.unlink(missing_ok=True)


def _dataset_target(project: Path | None):
    """The current project's loaded configuration, or None with an explanation."""

    from rnaseq.project import load_project

    if project is None:
        st.info("To use this data, choose a project in the sidebar or create one on **New Project**.")
        return None
    try:
        return load_project(project)
    except ProjectConfigError as exc:
        st.error(str(exc))
        return None


def datasets_page() -> None:
    from rnaseq.project import fastq_samplesheet_text, import_count_matrix

    config = gui_config()
    st.title("Datasets")
    st.caption("Browse the research data on this server. Only the data folders in the GUI settings are shown; FASTQ contents are never read here.")
    statuses = config.root_statuses()
    for item in statuses:
        if not item.available:
            st.error(f"`{item.path}`: {item.message}")
    roots = config.available_roots()
    if not roots:
        return
    data_root = config.project_parent.parent
    start = data_root / "fastq" if (data_root / "fastq").is_dir() and any((data_root / "fastq").iterdir()) else roots[0]
    folder = path_picker(config, "Dataset folder", "ds", want="dir", start=start)
    if folder is None:
        st.caption("Choose a folder: browse below and press *Use this folder*, or type its path.")
        return
    if st.button("Rescan", key="ds_rescan"):
        _scan.clear()
    try:
        scan = _scan(str(folder), tuple(map(str, roots)))
    except PathNotAllowed as exc:
        st.error(str(exc))
        return

    cols = st.columns(5)
    cols[0].metric("Files", len(scan.files))
    cols[1].metric("Total size", datasets.human_size(scan.total_bytes))
    cols[2].metric("FASTQ files", len(scan.of_kind("fastq")))
    cols[3].metric("Tables", sum(len(scan.of_kind(kind)) for kind in ("count_matrix", "metadata", "samplesheet", "contrasts", "table")))
    cols[4].metric("Unreadable", len(scan.unreadable))
    if scan.truncated:
        st.warning(f"Only the first {datasets.MAX_FILES} files are listed.")
    if scan.outside_links:
        st.warning("Symbolic links that point outside the allowed data folders are ignored: " + ", ".join(scan.outside_links[:10]))
    for problem in scan.unlisted:
        st.warning(f"Cannot list {problem}")
    if scan.unreadable:
        st.error("Files without read permission: " + ", ".join(item.relative for item in scan.unreadable[:10]))
    if scan.kind_summary():
        st.dataframe(pd.DataFrame(scan.kind_summary()), hide_index=True)
    with st.expander(f"All files ({len(scan.files)})", expanded=len(scan.files) <= 30):
        st.dataframe(pd.DataFrame([{"file": item.relative, "kind": datasets.KIND_LABELS[item.kind], "size": datasets.human_size(item.size),
                                    "readable": "yes" if item.readable else "NO", "link": "symlink" if item.symlink else ""} for item in scan.files[:2000]]),
                     hide_index=True, width="stretch")
    st.caption(f"Folders are listed {datasets.MAX_DEPTH} levels deep; run, work and output folders ({', '.join(sorted(datasets.SKIPPED_DIRS))}) are skipped.")

    project = current_project()
    tabs_wanted = [name for name, present in (("FASTQ", scan.fastq_folders), ("Count matrix", scan.of_kind("count_matrix")),
                                               ("Metadata", scan.of_kind("metadata")), ("Samplesheet", scan.of_kind("samplesheet"))) if present]
    if not tabs_wanted:
        st.info("No FASTQ files, count matrices, metadata or samplesheets were recognised in this folder.")
        return
    st.subheader("Recognise samples and use in a project")
    st.caption(f"Target project: **{project.name if project else 'none selected'}** (change it in the sidebar).")
    if project is not None:
        show_import_result(project)
    for tab, name in zip(st.tabs(tabs_wanted), tabs_wanted):
        with tab:
            if name == "FASTQ":
                _dataset_fastq(config, scan, project, fastq_samplesheet_text)
            elif name == "Count matrix":
                _dataset_counts(scan, project, import_count_matrix)
            elif name == "Metadata":
                _dataset_metadata(scan, project)
            else:
                _dataset_samplesheet(config, scan, project)


def _dataset_fastq(config: GuiConfig, scan: datasets.DatasetScan, project: Path | None, sheet_text) -> None:
    folders = scan.fastq_folders
    folder = folders[0] if len(folders) == 1 else Path(st.selectbox("FASTQ folder", [str(item) for item in folders], key="ds_fq_folder"))
    files = [item for item in scan.of_kind("fastq") if item.path.parent == folder]
    loaded = None
    if project is not None:
        try:
            from rnaseq.project import load_project
            loaded = load_project(project)
        except ProjectConfigError:
            loaded = None
    default = loaded.config.input.layout if loaded is not None and loaded.config.input.layout else datasets.guess_layout(files)
    layouts = [SequencingLayout.PAIRED_END, SequencingLayout.SINGLE_END]
    layout = st.radio("Layout", layouts, index=layouts.index(default), format_func=lambda value: value.value.replace("_", "-"), horizontal=True, key=f"ds_fq_layout::{folder}")
    summary, errors = datasets.recognise_fastq(folder, layout)
    st.caption(f"`{folder}` · {len(files)} FASTQ file(s), {datasets.human_size(sum(item.size or 0 for item in files))}. Samples follow the validator's naming rule `<sample>[_Lnnn]_R1|R2[_nnn].fastq.gz`.")
    show_issues(errors)
    if summary is None or not summary.records:
        return
    st.dataframe(pd.DataFrame([{"sample_id": record.sample_id, "lane": record.lane, "R1": record.fastq_1.name, "R2": record.fastq_2.name if record.fastq_2 else ""} for record in summary.records]), hide_index=True, width="stretch")
    st.caption(f"{len(summary.sample_ids)} sample(s), {len(summary.records)} FASTQ assignment(s).")
    loaded = _dataset_target(project)
    if loaded is None:
        return
    if loaded.config.input.type is not InputType.FASTQ:
        st.info(f"Project `{project.name}` takes a raw-count matrix, not FASTQ files.")
        return
    if loaded.config.input.layout is not layout:
        st.warning(f"Project `{project.name}` is {loaded.config.input.layout.value}; recognise the files as that layout to import them.")
        return
    if errors:
        st.warning("Fix the FASTQ naming problems above before importing.")
        return
    plan = datasets.fastq_copy_plan(summary.records, loaded.input_path)
    copy_plan_notice(plan)
    sheet = sheet_text(summary.records, loaded.config.upstream.strandedness or "auto")
    if st.button(f"Copy {len(summary.records)} FASTQ assignment(s) into {project.name}", key="ds_fq_import", disabled=not plan.fits, type="primary"):
        import_fastq(config, project, summary.records, sheet_text=sheet)


def _dataset_counts(scan: datasets.DatasetScan, project: Path | None, importer) -> None:
    files = scan.of_kind("count_matrix")
    chosen = files[0] if len(files) == 1 else files[[item.relative for item in files].index(st.selectbox("Count matrix", [item.relative for item in files], key="ds_counts_file"))]
    summary, issues = _recognise_counts(str(chosen.path), chosen.path.stat().st_mtime)
    if summary is not None:
        st.success(f"`{chosen.relative}`: {summary.gene_count} genes × {len(summary.sample_ids)} samples: {', '.join(summary.sample_ids[:20])}{' …' if len(summary.sample_ids) > 20 else ''}")
    show_issues(issues)
    loaded = _dataset_target(project)
    if loaded is None or summary is None:
        return
    if loaded.config.input.type is not InputType.RAW_COUNTS:
        st.info(f"Project `{project.name}` takes FASTQ files, not a count matrix.")
        return
    st.caption(f"Importing copies this file ({datasets.human_size(chosen.size)}) to `{loaded.input_path}`.")
    replace = loaded.input_path.exists() and st.checkbox("Replace the project's existing count matrix", key="ds_counts_replace")
    if st.button(f"Use this count matrix in {project.name}", key="ds_counts_import", type="primary", disabled=loaded.input_path.exists() and not replace):
        try:
            importer(project, chosen.path, replace=bool(replace))
        except ProjectCreationError as exc:
            st.error(str(exc))
        else:
            st.success(f"Copied to `{loaded.input_path}`. Next: describe the samples on Input & Metadata.")
            if "inputs" in PAGES:
                st.page_link(PAGES["inputs"], label="Input & Metadata →")


def _dataset_metadata(scan: datasets.DatasetScan, project: Path | None) -> None:
    files = scan.of_kind("metadata")
    chosen = files[0] if len(files) == 1 else files[[item.relative for item in files].index(st.selectbox("Metadata file", [item.relative for item in files], key="ds_meta_file"))]
    header, rows = datasets.table_preview(chosen.path, 50)
    st.dataframe(pd.DataFrame(rows, columns=header), hide_index=True, width="stretch")
    loaded = _dataset_target(project)
    if loaded is None:
        return
    if loaded.config.project.preset is Preset.QC:
        st.info("QC-only projects do not use metadata.")
        return
    report = validated(project)
    if report.config is None:
        show_issues(report.errors)
        return
    parsed_header, parsed_rows = services.parse_table_text(chosen.path.read_text(encoding="utf-8-sig"))
    columns, rows = _table_rows(editor_table(project, report))
    change = services.merge_metadata_table(columns, rows, parsed_header, parsed_rows)
    st.caption(f"Merged by sample_id into the metadata being edited for {project.name} (unsaved edits included).")
    ok = metadata_change_preview(change, report, f"ds_meta_confirm::{project}::{chosen.relative}")
    if st.button(f"Load into the metadata editor of {project.name}", key="ds_meta_load", disabled=not ok):
        _apply_metadata_change(project, change, f"{chosen.relative} applied")
        st.session_state["_ds_meta_loaded"] = _change_summary(change)
        st.rerun()
    loaded_summary = st.session_state.pop("_ds_meta_loaded", None)
    if loaded_summary:
        st.success(f"Loaded ({loaded_summary}), not saved yet. Check the table and press *Save metadata* on Input & Metadata.")
        if "inputs" in PAGES:
            st.page_link(PAGES["inputs"], label="Input & Metadata →")


def _dataset_samplesheet(config: GuiConfig, scan: datasets.DatasetScan, project: Path | None) -> None:
    files = scan.of_kind("samplesheet")
    chosen = files[0] if len(files) == 1 else files[[item.relative for item in files].index(st.selectbox("Samplesheet", [item.relative for item in files], key="ds_sheet_file"))]
    view = datasets.recognise_samplesheet(config, chosen.path)
    if view.error:
        st.error(view.error)
        return
    st.dataframe(pd.DataFrame([{"sample": record.sample, "fastq_1": str(record.fastq_1), "fastq_2": str(record.fastq_2 or ""), "strandedness": record.strandedness} for record in view.records]), hide_index=True, width="stretch")
    st.caption(f"{len(view.records)} row(s), {view.layout.value if view.layout else '?'}, strandedness {view.strandedness}.")
    if view.outside:
        st.error("Rows refer to FASTQ files outside the allowed data folders; this samplesheet cannot be imported here: " + ", ".join(view.outside[:5]))
        return
    loaded = _dataset_target(project)
    if loaded is None:
        return
    if loaded.config.input.type is not InputType.FASTQ:
        st.info(f"Project `{project.name}` takes a raw-count matrix, not FASTQ files.")
        return
    plan = datasets.fastq_copy_plan(view.records, loaded.input_path)
    copy_plan_notice(plan)
    if st.button(f"Import this samplesheet into {project.name}", key="ds_sheet_import", disabled=not plan.fits, type="primary"):
        import_fastq(config, project, view.records, sheet_path=chosen.path)


# --------------------------------------------------------------------- B. New Project

OTHER_REFERENCE = "Other (advanced)"


REFERENCE_REGISTRY_HELP = (
    "Managed references are offered from this workstation's registry (`{registry}`). "
    "A prepared reference folder (with `reference_manifest.yaml`) is added with `rnaseq reference register /absolute/reference-root`; "
    "an index is built with `rnaseq reference prepare-hisat2 ROOT` (HISAT2) or `rnaseq reference prepare ROOT` (Salmon). "
    "See docs/runtime.md. Nothing is downloaded or built from the GUI."
)


def reference_not_prepared(species: str, method: str, available: services.ReferenceAvailability) -> None:
    """Say that no registered reference fits, and why, instead of offering an empty or unrelated path."""

    from rnaseq.references import reference_registry_path

    backend = services.BACKEND_LABELS.get(method, method)
    lines = [f"**{species} {backend} reference not prepared** — no registered reference of this species has a usable {backend} index."]
    for item, reason in available.unusable:
        lines.append(f"- `{item.root}` ({item.provider} {item.release} / {item.assembly_identity}): {reason}")
    st.warning("\n".join(lines))
    st.caption(REFERENCE_REGISTRY_HELP.format(registry=reference_registry_path()))


def _mtime(path: Path) -> float | None:
    try:
        return path.stat().st_mtime
    except OSError:
        return None


@st.cache_data(max_entries=16, show_spinner="Checking the reference manifest ...")
def _cached_reference_check(root: str, manifest: str, species: str, method: str, mtime: float | None) -> str | None:
    return services.check_project_reference({"source": "local", "root": root, "manifest": manifest}, species, method)


def _other_reference(species: str, method: str) -> dict[str, object] | None:
    """The fields of an unregistered reference, shown as soon as *Other* is chosen.

    Returns None (with a note, not an error) until the required fields are filled.
    """

    box = st.container(border=True)
    box.markdown("**Other reference**")
    source = box.radio("Reference type", ["local", "custom"], horizontal=True, key="np_ref_source",
                       format_func=lambda value: {"local": "Managed reference folder (with a reference manifest)", "custom": "Genome FASTA + annotation GTF"}[value])
    if source == "local":
        root = box.text_input("Managed reference folder", key="np_ref_root", help="The folder that holds the reference manifest, as prepared by `rnaseq reference`.")
        manifest = box.text_input("Reference manifest file name", value="reference_manifest.yaml", key="np_ref_manifest")
        if not root.strip():
            box.caption("Enter the managed reference folder to use this reference.")
            return None
        reference = wizard.reference_options(source="local", species=wizard.species_value(species), method=method, local_root=Path(root.strip()),
                                             local_manifest=manifest.strip() or None, fasta=None, gtf=None, transcript_fasta=None, salmon_index=None, hisat2_index=None)
        problem = _cached_reference_check(str(reference["root"]), str(reference["manifest"]), wizard.species_value(species).value, method, _mtime(Path(str(reference["root"])) / str(reference["manifest"])))
        if problem:
            box.error(problem)
            return None
        return reference
    fasta = box.text_input("Genome FASTA", key="np_fasta").strip()
    gtf = box.text_input("Annotation GTF", key="np_gtf").strip()
    hisat2 = box.text_input("Prepared HISAT2 index (optional)", key="np_hisat2").strip()
    missing = [name for name, value in (("genome FASTA", fasta), ("annotation GTF", gtf)) if not value]
    if missing:
        box.caption(f"Enter the {' and the '.join(missing)} to use this reference.")
        return None
    return wizard.reference_options(source="custom", species=wizard.species_value(species), method=method, local_root=None, local_manifest=None,
                                    fasta=Path(fasta), gtf=Path(gtf), transcript_fasta=None, salmon_index=None, hisat2_index=Path(hisat2) if hisat2 else None)


def new_project_page() -> None:
    from rnaseq.execution import detect_local_resource_capacity, validate_local_execution_budget
    from rnaseq.project import create_project
    from rnaseq.resource_policy import AUTO

    config = gui_config()
    st.title("New project")
    st.caption("Same choices as the interactive `rnaseq new` wizard. The project is created as a template (scaffold); data and metadata are added on the next page.")

    name = st.text_input("Project name", key="np_name", help="Letters, digits, '.', '_' or '-'; also the folder name.")
    parent_text = st.text_input("Location (parent folder)", value=str(config.project_parent), key="np_parent")
    species = st.selectbox("Species", [value for value, _ in wizard.SPECIES_CHOICES], format_func=dict(wizard.SPECIES_CHOICES).get, key="np_species")
    input_type = st.radio("Input", [value for value, _ in wizard.INPUT_TYPE_CHOICES], format_func=lambda value: f"{value} — {dict(wizard.INPUT_TYPE_CHOICES)[value]}", key="np_input", horizontal=True)

    layout = preprocessing = method = strandedness = None
    reference: dict[str, object] | None = None
    if input_type == "fastq":
        cols = st.columns(3)
        methods = [value for value, _ in wizard.METHOD_CHOICES]
        method = cols[0].selectbox("Quantification", methods, index=methods.index(workspace.DEFAULT_METHOD), format_func=lambda value: workspace.METHOD_LABELS.get(value, value), key="np_method")
        layout = cols[1].selectbox("Sequencing layout", [value for value, _ in wizard.LAYOUT_CHOICES], format_func=lambda value: f"{value} ({dict(wizard.LAYOUT_CHOICES)[value]})", key="np_layout")
        strand_options = [value for value, _ in wizard.strandedness_choices(method)]
        if method == "hisat2_featurecounts":
            strandedness = cols[2].selectbox("Strandedness", strand_options, index=None, placeholder="Choose explicitly", key="np_strand_hisat2")
        else:
            strandedness = cols[2].selectbox("Strandedness", strand_options, index=strand_options.index(wizard.default_strandedness(method)), format_func=lambda value: f"{value} ({dict(wizard.strandedness_choices(method))[value]})", key="np_strand_salmon")
        st.caption(BACKEND_HELP.get(method, ""))
        if method == "hisat2_featurecounts":
            with st.expander("Which strandedness? (never guessed for HISAT2 + featureCounts)", expanded=strandedness is None):
                st.markdown(STRANDEDNESS_HELP)
        preprocessing = st.radio("FASTQ preprocessing", [value for value, _ in wizard.PREPROCESSING_CHOICES], format_func=lambda value: f"{value} — {dict(wizard.PREPROCESSING_CHOICES)[value]}",
                                 key="np_pre", horizontal=True,
                                 help="raw: nf-core/rnaseq trims adapters and low-quality bases before quantification. "
                                      "pretrimmed: the reads were already trimmed (for example by the sequencing provider) and are used as supplied.")
        species_name = wizard.species_value(species).value
        try:
            available = services.reference_availability(species_name, method)
        except Exception as exc:  # a broken registry must not block the other reference routes
            available = services.ReferenceAvailability([], [], str(exc))
        if available.error:
            st.warning(f"Reference registry could not be read: {available.error}")
        registered = available.compatible
        if not registered:
            reference_not_prepared(species_name, method, available)
        labels = {f"{ref.provider} {ref.release} / {ref.assembly_identity} ({ref.root})": ref for ref in registered}
        other = ["iGenomes (Salmon only)"] if method == "salmon" else []
        # Keyed by species and backend: a choice made for one combination (for
        # example "Other", the only option while Human has no registered
        # reference) is never carried over to another, where it would hide the
        # registered reference that applies there.
        choice = st.selectbox("Reference", [*labels, *other, OTHER_REFERENCE], key=f"np_reference::{species}::{method}",
                              help="Only references registered with `rnaseq reference register` that have a usable index for this species and quantification are listed.")
        try:
            if choice in labels:
                root, manifest = wizard.managed_reference_settings(labels[choice])
                reference = wizard.reference_options(source="local", species=wizard.species_value(species), method=method, local_root=root, local_manifest=manifest, fasta=None, gtf=None, transcript_fasta=None, salmon_index=None, hisat2_index=None)
            elif choice in other:
                reference = wizard.reference_options(source="igenomes", species=wizard.species_value(species), method=method, local_root=None, local_manifest=None, fasta=None, gtf=None, transcript_fasta=None, salmon_index=None, hisat2_index=None)
            else:
                reference = _other_reference(species, method)
        except ProjectCreationError as exc:
            st.error(str(exc))

    preset = st.radio("Analysis scope", [value for value, _ in wizard.preset_choices(input_type)], format_func=lambda value: f"{value} — {dict(wizard.preset_choices(input_type))[value]}", key=f"np_preset_{input_type}", index=1 if input_type == "fastq" else 0)
    design_type, pair_id, batch, enrichment, gene_id_type = "two_group", None, False, None, None
    if preset != "qc":
        design_type = st.selectbox("Experimental design", [value for value, _ in wizard.DESIGN_CHOICES], format_func=lambda value: f"{value} — {dict(wizard.DESIGN_CHOICES)[value]}", key="np_design")
        st.caption(DESIGN_HELP.get(design_type, ""))
        if design_type == "paired_two_group":
            pair_id = st.text_input("Pairing metadata column", value=wizard.DEFAULT_PAIR_ID, key="np_pair", help="The metadata column naming the biological pairing unit (patient, donor, cell line ...).")
        batch = st.checkbox("There is a known technical batch (adds a 'batch' column; the model adjusts for it)", key="np_batch")
        if design_type != "paired_two_group" or (pair_id or "").strip():
            from rnaseq.project import derive_design

            formula, _types, design_columns = derive_design(DesignType(design_type), (pair_id or "").strip() or None, batch)
            st.markdown(f"DESeq2 design formula: `{formula}` · metadata columns: `sample_id`, " + ", ".join(f"`{column}`" for column in design_columns))
        if preset == "L2":
            enrichment = st.checkbox(f"Run functional enrichment: {describe_enrichment(DEFAULT_L2_ENRICHMENT)}", value=True, key="np_enrichment")
        if preset == "L1" or enrichment:
            gene_id_type = st.selectbox("Gene identifier type in the count matrix", [value for value, _ in wizard.GENE_ID_TYPE_CHOICES], format_func=lambda value: f"{value} — {dict(wizard.GENE_ID_TYPE_CHOICES)[value]}", key="np_geneid")

    suggestion = wizard.resource_suggestion()
    with st.expander("Local execution resources"):
        st.caption(f"Detected {suggestion.logical_cpus or '?'} CPUs / {suggestion.memory_gib or '?'} GiB; `auto` gives {suggestion.auto_cpus} CPUs / {suggestion.auto_memory_gib} GiB on this machine (resolved again on each run).")
        cpus = st.number_input("CPU limit (0 = auto)", min_value=0, value=0, step=1, key="np_cpus")
        memory = st.number_input("Memory limit in GiB (0 = auto)", min_value=0, value=0, step=1, key="np_memory")

    missing = [item for item, absent in (
        ("a project name", not name.strip()),
        ("the strandedness (HISAT2 + featureCounts)", input_type == "fastq" and method == "hisat2_featurecounts" and not strandedness),
        ("a reference", input_type == "fastq" and reference is None),
        ("the pairing column", preset != "qc" and design_type == "paired_two_group" and not (pair_id or "").strip()),
    ) if absent]
    if missing:
        st.caption("Still needed before *Create project*: " + "; ".join(missing) + ".")
    else:
        st.caption(f"*Create project* writes a new folder `{parent_text.rstrip('/')}/{name.strip()}` (project.yaml and empty metadata/contrast templates). "
                   "Nothing else is changed; data and metadata are added on the next page.")
    if st.button("Create project", type="primary", key="np_create"):
        try:
            parent = config.resolve_allowed(parent_text)
            if not name:
                raise ProjectCreationError("Enter a project name.")
            if input_type == "fastq" and method == "hisat2_featurecounts" and not strandedness:
                raise ProjectCreationError("Choose the strandedness explicitly for HISAT2 + featureCounts.")
            if input_type == "fastq" and reference is None:
                raise ProjectCreationError("Choose a reference, or fill in the required fields of the other reference.")
            if input_type == "fastq":
                # Checked again here (the selection may predate a species or backend change); the validator checks it once more at run time.
                problem = services.check_project_reference(reference, wizard.species_value(species).value, str(method))
                if problem:
                    raise ProjectCreationError(problem)
            selected_cpus, selected_memory = int(cpus) or AUTO, int(memory) or AUTO
            validate_local_execution_budget(selected_cpus, selected_memory, detect_local_resource_capacity())
            raw = input_type == "raw_counts"
            target = create_project(
                project_name=name, destination=parent, species=wizard.species_value(species), preset=Preset.QC if preset == "qc" else Preset(preset),
                design_type=DesignType(design_type), input_type=InputType(input_type),
                layout=None if raw else SequencingLayout(layout), preprocessing=FastqPreprocessing(preprocessing or "raw"),
                strandedness="auto" if raw else str(strandedness), quantification_method="salmon" if raw else str(method),
                reference={"source": "igenomes", "genome": None} if raw else reference, scaffold=True, pair_id=pair_id,
                execution={"profile": "local", "max_cpus": selected_cpus, "max_memory_gb": selected_memory}, batch=batch,
                enrichment=None if enrichment is None else (DEFAULT_L2_ENRICHMENT if enrichment else ()), gene_id_type=gene_id_type,
            )
        except (PathNotAllowed, ProjectCreationError, ExecutionPreflightError, ValueError) as exc:
            st.error(str(exc))
            return
        except OSError as exc:
            st.error(f"System error: {exc}")
            return
        select_project(target)
        st.success(f"Created project `{target}`.")
        if "inputs" in PAGES:  # absent when a page runs on its own (AppTest)
            st.page_link(PAGES["datasets"], label="Next: choose its data on Datasets →")
            st.page_link(PAGES["inputs"], label="or add data and metadata on Input & Metadata →")


# --------------------------------------------------------------------- C. Input & Metadata

def _input_section(config: GuiConfig, project: Path, report) -> None:
    from rnaseq.project import fastq_samplesheet_text, import_count_matrix

    assert report.config is not None and report.loaded is not None
    st.subheader("1. Input data")
    st.caption("Tip: the **Datasets** page lists everything in the data folders and can assign it to this project.")
    show_import_result(project)
    if report.config.input.type is InputType.RAW_COUNTS:
        if report.counts is not None:
            st.success(f"Count matrix `{report.loaded.input_path}`: {report.counts.gene_count} genes × {len(report.counts.sample_ids)} samples.")
        else:
            st.warning(f"No usable count matrix at `{report.loaded.input_path}` yet.")
        source = path_picker(config, "Count matrix file (CSV, gene_id first column, integer counts)", "counts", want="file", suffixes=(".csv", ".tsv", ".txt"))
        replace = report.loaded.input_path.exists() and st.checkbox("Replace the existing count matrix", key="counts_replace")
        if source is not None and st.button("Import count matrix", key="counts_import"):
            try:
                import_count_matrix(project, source, replace=bool(replace))
                st.success(f"Copied {source.name} to {report.loaded.input_path}.")
                st.rerun()
            except ProjectCreationError as exc:
                st.error(str(exc))
        return

    layout = report.config.input.layout or SequencingLayout.PAIRED_END
    if report.fastq is not None and report.fastq.records:
        st.success(f"FASTQ input `{report.loaded.input_path}`: {len(report.fastq.sample_ids)} samples, {len(report.fastq.records)} FASTQ assignments ({layout.value}).")
        origin = workspace.input_origin(project)
        if origin:
            st.caption("Copied from: " + ", ".join(f"`{item}`" for item in origin))
    else:
        st.warning(f"No FASTQ files in `{report.loaded.input_path}` yet.")
    mode = st.radio("Add FASTQ files", ["Detect in a folder", "Import an existing samplesheet"], horizontal=True, key="fq_mode",
                    help="Detect: choose a folder; sample IDs and R1/R2 are recognised from the file names. "
                         "Samplesheet: a CSV with sample,fastq_1,fastq_2,strandedness. Both show a preview first.")
    st.caption("Choosing a folder or samplesheet only **previews** the assignments; nothing is written until you press the *Copy* / *Import* button, "
               "which copies the FASTQ files into the project (the originals are not changed).")
    if mode == "Detect in a folder":
        folder = path_picker(config, "FASTQ folder", "fastq_dir", want="dir", start=config.project_parent.parent / "fastq")
        if folder is None:
            return
        summary, errors = services.discover_fastq(folder, layout)
        show_issues(errors)
        if summary is None or not summary.records:
            return
        st.markdown(f"**Preview** — {len(summary.sample_ids)} sample(s), {len(summary.records)} FASTQ assignment(s) detected in `{folder}` (not yet copied):")
        st.dataframe(pd.DataFrame([{"sample_id": record.sample_id, "lane": record.lane, "R1": record.fastq_1.name, "R2": record.fastq_2.name if record.fastq_2 else ""} for record in summary.records]), hide_index=True, width="stretch")
        sheet = fastq_samplesheet_text(summary.records, report.config.upstream.strandedness or "auto")
        st.code(sheet, language=None)
        st.download_button("Download samplesheet.csv", sheet, file_name="samplesheet.csv", key="fq_download")
        if errors:
            return
        plan = datasets.fastq_copy_plan(summary.records, report.loaded.input_path)
        copy_plan_notice(plan)
        if st.button(f"Copy these {len(summary.records)} FASTQ assignment(s) into the project", key="fq_import", disabled=not plan.fits):
            # The import keeps its own copy as planning/imported_fastq_samplesheet.csv.
            import_fastq(config, project, summary.records, sheet_text=sheet)
    else:
        sheet_path = path_picker(config, "Samplesheet CSV (sample,fastq_1,fastq_2,strandedness)", "fq_sheet", want="file", suffixes=(".csv",))
        if sheet_path is None:
            return
        view = datasets.recognise_samplesheet(config, sheet_path)
        if view.error:
            st.error(view.error)
            return
        if view.outside:
            st.error("Rows refer to FASTQ files outside the allowed data folders: " + ", ".join(view.outside[:5]))
            return
        plan = datasets.fastq_copy_plan(view.records, report.loaded.input_path)
        copy_plan_notice(plan)
        if st.button("Import samplesheet (copies the FASTQ files)", key="fq_sheet_import", disabled=not plan.fits):
            import_fastq(config, project, view.records, sheet_path=sheet_path)


IMPORT_RESULT_KEY = "_fastq_import_result"


def _import_samplesheet(project: Path, sheet: Path, importer, plan: datasets.CopyPlan | None = None) -> None:
    """Run the existing importer; it reports no per-file progress, so the page says it is blocked and then what was copied."""

    size = f" ({plan.files} file(s), {datasets.human_size(plan.bytes)})" if plan is not None else ""
    notice = st.warning(f"Copying FASTQ files{size} into the project. This blocks the GUI until every file is copied; "
                        "keep this tab open and do not start another import. The importer does not report per-file progress.")
    started = datetime.now()
    try:
        with st.spinner("Copying FASTQ files into the project ..."):
            records = importer(project, sheet)
    except (ProjectCreationError, ProjectConfigError) as exc:
        notice.empty()
        st.error(f"Not imported: {exc}")
        return
    except OSError as exc:
        notice.empty()
        st.error(f"System error while copying: {exc}. Files already copied stay in the project input folder; check it before importing again.")
        return
    seconds = int((datetime.now() - started).total_seconds())
    st.session_state[IMPORT_RESULT_KEY] = (str(project), f"Copied the FASTQ files of {len(records)} samplesheet row(s){size} into the project in {seconds} s. "
                                           "The originals were not changed; the samplesheet is kept as planning/imported_fastq_samplesheet.csv.")
    st.rerun()


def show_import_result(project: Path) -> None:
    result = st.session_state.get(IMPORT_RESULT_KEY)
    if result and result[0] == str(project):
        st.success(result[1])
        st.session_state.pop(IMPORT_RESULT_KEY, None)


def _metadata_state_key(project: Path) -> str:
    return f"metadata::{project}"


def _set_editor_table(project: Path, frame: pd.DataFrame) -> None:
    """Make ``frame`` the table being edited, in a fresh editor (no pending edits are applied on top of it)."""

    key = _metadata_state_key(project)
    frame = frame.reset_index(drop=True)
    st.session_state[key] = frame  # the data handed to st.data_editor
    st.session_state[f"{key}::working"] = frame  # what the user sees: the data plus the editor's own edits
    st.session_state[f"{key}::version"] = st.session_state.get(f"{key}::version", 0) + 1
    st.session_state.pop(f"{key}::result", None)  # a save message would describe an older table


def editor_table(project: Path, report) -> pd.DataFrame:
    """The metadata table being edited, unsaved edits included.

    Initialised from metadata.csv plus the input's sample IDs. The editor's
    cell edits are copied to ``<key>::working`` on every rerun, so they survive
    bulk assignment, pasting, importing and visits to other pages.
    """

    key = _metadata_state_key(project)
    if f"{key}::working" in st.session_state:
        return st.session_state[f"{key}::working"]
    columns = services.metadata_columns(report)
    rows = services.merge_sample_ids(services.read_metadata(report), columns, services.input_sample_ids(report))
    frame = pd.DataFrame(rows, columns=columns)
    _set_editor_table(project, frame)
    return frame


def _table_rows(table: pd.DataFrame) -> tuple[list[str], list[dict]]:
    return [str(column) for column in table.columns], table.to_dict("records")


def _apply_metadata_change(project: Path, change: services.MetadataChange, what: str) -> None:
    _set_editor_table(project, pd.DataFrame(change.rows, columns=change.columns))
    st.session_state[f"{_metadata_state_key(project)}::notice"] = f"{what}: {_change_summary(change)}. Not saved yet."


def _change_summary(change: services.MetadataChange) -> str:
    parts = []
    if change.filled:
        parts.append(f"filled {change.filled} empty value(s)")
    if change.overwritten:
        parts.append(f"replaced {len(change.overwritten)} value(s)")
    if change.added:
        parts.append(f"added {len(change.added)} sample row(s) ({', '.join(change.added[:8])}{' …' if len(change.added) > 8 else ''})")
    if change.new_columns:
        parts.append("added column(s) " + ", ".join(change.new_columns))
    return ", ".join(parts) or "no change"


def metadata_change_preview(change: services.MetadataChange, report, confirm_key: str, *, merged: bool = True) -> bool:
    """Show what applying ``change`` does to the table being edited; True when it may be applied.

    Replacing existing non-empty values needs an explicit tick.
    """

    for error in change.errors:
        st.error(error)
    if change.errors:
        return False
    if not change.changes:
        st.caption("Nothing would change: the table already has these values.")
        return False
    text = _change_summary(change).replace("filled", "fills").replace("replaced", "replaces").replace("added", "adds")
    note = " Samples and columns that are not in the new values keep their current values; empty cells in the new values change nothing." if merged else ""
    st.info(f"Applying {text}.{note}")
    known = services.input_sample_ids(report)
    unknown = [sample for sample in change.added if known and sample not in known]
    if unknown:
        st.warning(f"{len(unknown)} added sample(s) are not samples of the input (count matrix or FASTQ): {', '.join(unknown[:10])}. Saving will report the disagreement.")
    if change.overwritten:
        st.warning(f"{len(change.overwritten)} existing value(s) would be replaced:")
        st.dataframe(pd.DataFrame([{"sample_id": sample, "column": column, "current value": old, "new value": new} for sample, column, old, new in change.overwritten]),
                     hide_index=True, width="stretch")
        return st.checkbox(f"Replace these {len(change.overwritten)} existing value(s)", key=confirm_key)
    return True


def _metadata_section(project: Path, report) -> None:
    assert report.config is not None
    st.subheader("2. Sample metadata")
    if report.config.project.preset is Preset.QC:
        st.info("QC-only projects make no statistical design claim; metadata is not used.")
        return
    design_vars = services.design_variables(report)
    st.caption(
        f"Design formula `{report.config.design.formula}`; required columns: sample_id, {', '.join(design_vars)}. "
        "Conditions are entered by you; they are never derived from file names."
    )
    key = _metadata_state_key(project)
    current = editor_table(project, report)
    editor_key = f"{key}::editor::{st.session_state.get(f'{key}::version', 0)}"
    if editor_key not in st.session_state and not current.equals(st.session_state[key]):
        # Streamlit dropped the editor's own state (the page was left); its
        # edits, kept in ::working, become the editor's data.
        st.session_state[key] = current.reset_index(drop=True)
    # Always drawn (empty or not) so the tabs below keep their position and the
    # chosen tab stays open after an Apply.
    notice_box = st.empty()
    notice = st.session_state.pop(f"{key}::notice", None)
    if notice:
        notice_box.success(notice)

    tools_box = st.container()  # drawn above the table, filled after it so the tools act on the current edits
    base: pd.DataFrame = st.session_state[key]
    edited = st.data_editor(
        base, num_rows="dynamic", hide_index=True, width="stretch", key=editor_key,
        column_config={column: st.column_config.TextColumn(column, required=column in ("sample_id", *design_vars)) for column in base.columns},
    )
    st.session_state[f"{key}::working"] = edited
    columns, rows = _table_rows(edited)
    difference = services.metadata_difference(services.metadata_file_columns(report), services.read_metadata(report), columns, rows)
    if difference.unsaved:
        st.warning(f"**Unsaved changes** — {difference.describe()}. Press *Save metadata* to write them to metadata.csv.")
    else:
        st.caption("No unsaved changes: the table matches metadata.csv.")

    with tools_box:
        tools = st.tabs(["Paste table (TSV / Excel)", "Import CSV", "Assign values to samples", "Reload"])
        with tools[0]:
            pasted = st.text_area("Paste rows copied from Excel (tab-separated) or CSV, header row first (must include sample_id)", key=f"{key}::paste", height=150)
            if pasted.strip():
                header, new_rows = services.parse_table_text(pasted)
                change = services.merge_metadata_table(columns, rows, header, new_rows)
                ok = metadata_change_preview(change, report, f"{key}::paste_confirm")
                if st.button("Apply pasted table", key=f"{key}::paste_apply", disabled=not ok):
                    _apply_metadata_change(project, change, "Pasted table applied")
                    st.rerun()
        with tools[1]:
            uploaded = st.file_uploader("Metadata CSV/TSV (must include sample_id)", type=["csv", "tsv", "txt"], key=f"{key}::upload")
            if uploaded is not None:
                header, new_rows = services.parse_table_text(uploaded.getvalue().decode("utf-8-sig"))
                change = services.merge_metadata_table(columns, rows, header, new_rows)
                ok = metadata_change_preview(change, report, f"{key}::upload_confirm")
                if st.button("Load into the editor", key=f"{key}::upload_apply", disabled=not ok):
                    _apply_metadata_change(project, change, f"{uploaded.name} applied")
                    st.rerun()
        with tools[2]:
            # A fresh widget key after each Apply clears the previous selection.
            bulk = f"{key}::bulk::{st.session_state.get(f'{key}::version', 0)}"
            cols = st.columns([3, 2, 2, 1])
            samples = cols[0].multiselect("Samples", [str(value) for value in edited["sample_id"].tolist() if _filled(value)], key=f"{bulk}::samples")
            column = cols[1].selectbox("Column", [column for column in columns if column != "sample_id"], key=f"{bulk}::column")
            value = cols[2].text_input("Value", key=f"{bulk}::value")
            change = services.assign_metadata_value(columns, rows, samples, column, value) if samples and column else None
            ok = change is not None and metadata_change_preview(change, report, f"{bulk}::confirm", merged=False)
            if cols[3].button("Apply", key=f"{bulk}::apply", disabled=not ok, help=None if samples and column else "Choose samples and a column first."):
                _apply_metadata_change(project, change, f"{column} = {value!r} applied")
                st.rerun()
        with tools[3]:
            confirm = not difference.unsaved or st.checkbox(f"Discard the unsaved changes ({difference.describe()})", key=f"{key}::reload_confirm")
            if st.button("Reload metadata.csv", key=f"{key}::reload", disabled=not confirm):
                fresh = services.merge_sample_ids(services.read_metadata(report), services.metadata_columns(report), services.input_sample_ids(report))
                _set_editor_table(project, pd.DataFrame(fresh, columns=services.metadata_columns(report)))
                st.rerun()

    if report.config.design.type is DesignType.PAIRED_TWO_GROUP and report.config.design.pair_id:
        pairing_check(columns, rows, report)

    if st.button("Save metadata", type="primary", key=f"{key}::save", help="Writes the table to metadata.csv; it is validated first and not written if it has errors."):
        result = services.save_metadata(project, columns, [{column: services.clean_cell(row.get(column)) for column in columns} for row in rows])
        if result.saved:
            _set_editor_table(project, edited)
        st.session_state[f"{key}::result"] = (result, columns)
        st.rerun()
    stored = st.session_state.get(f"{key}::result")
    if stored is not None:
        result, saved_columns = stored
        if result.saved:
            st.success(f"Saved {report.loaded.metadata_path}.")
        else:
            st.error("metadata.csv was not changed; fix these first:")
        by_column: dict[str, list] = {}
        for issue in [*result.errors, *result.incomplete, *result.agreement]:
            by_column.setdefault(services.issue_column(issue, saved_columns), []).append(issue)
        for column, issues in by_column.items():
            st.markdown(f"**{column}**")
            show_issues(issues)
        if result.saved and not (result.incomplete or result.agreement):
            st.caption("All design values are filled. Next: Contrasts (L2) or Run.")


def pairing_check(columns: list[str], rows: list[dict], report) -> None:
    """Live check of the pairs in the table being edited (saved or not), by the validator's rule for paired designs."""

    from rnaseq.validators import condition_factor

    pair_column = report.config.design.pair_id
    factor = condition_factor(report.formula_variables, report.design_variable_types) or "condition"
    preview = services.pairing_preview(columns, rows, pair_column, factor)
    if preview is None:
        st.warning(f"Pairing check: the table needs the columns `{pair_column}` and `{factor}`.")
        return
    filled = preview.pairs or preview.blank_pair
    if not filled:
        st.caption(f"Pairing check: fill `{pair_column}` and `{factor}` to see the pairs.")
        return
    headline = (f"Pairing check (`{pair_column}` × `{factor}`, table as shown, saved or not): {preview.complete} complete pair(s) of "
                f"{len(preview.pairs)}" + (f"; {len(preview.blank_pair)} sample(s) without a `{pair_column}` value" if preview.blank_pair else ""))
    (st.success if preview.problems == 0 and len(preview.levels) == 2 else st.warning)(headline + ".")
    with st.expander("Pairs", expanded=preview.problems > 0):
        st.dataframe(pd.DataFrame([{pair_column: row.pair, **{level: ", ".join(row.samples.get(level, ())) for level in preview.levels}, "status": row.status}
                                   for row in preview.pairs]), hide_index=True, width="stretch")
        if preview.blank_pair:
            st.caption(f"Without a `{pair_column}` value: " + ", ".join(preview.blank_pair))
        st.caption("Each pair needs exactly one sample of each of the two condition levels (the rule `rnaseq validate` applies; "
                   "L2 needs at least two complete pairs). The validator remains the authority when you save and run.")


def _filled(value) -> bool:
    return bool(services.clean_cell(value))


def inputs_page() -> None:
    config = gui_config()
    st.title("Input & metadata")
    project = require_project()
    if project is None:
        return
    report = validated(project)
    if report.config is None:
        show_issues(report.errors)
        return
    st.caption(f"`{project}` · {services.scope_label(report)} · input {report.config.input.type.value}")
    workflow_bar(project, report, "inputs")
    _input_section(config, project, report)
    _metadata_section(project, report)


# --------------------------------------------------------------------- D. Contrast designer

def contrasts_page() -> None:
    from rnaseq.design_setup import (
        DesignRejected, DesignSetupError, apply_design, build_proposal, design_blockers, design_enrichment_selection,
        design_needs_gene_id_type, pin_kegg_snapshot, residual_df_blocker, resolve_design_enrichment,
    )

    gui_config()
    st.title("Contrasts")
    project = require_project()
    if project is None:
        return
    report = validated(project)
    config = report.config
    if config is None:
        show_issues(report.errors)
        return
    workflow_bar(project, report, "contrasts")
    if config.project.preset is Preset.QC:
        st.info("This project's scope is QC-only; there is no statistical design to configure.")
        return
    view = services.design_view(report)
    st.subheader("Design")
    st.markdown(
        f"- Design type: **{config.design.type.value}**\n- DESeq2 formula: `{view.formula}`\n"
        + (f"- Pairing column: `{view.pair_id}` (each comparison is made within pairs)\n" if view.paired else "")
        + f"- Comparison factor: `{view.factor}`"
    )
    blockers = design_blockers(report)
    if view.levels:
        st.dataframe(pd.DataFrame([{"level": level, "samples": view.counts[level]} for level in view.levels]), hide_index=True)
    if blockers:
        st.error("The metadata/design is not valid yet; fix these on the Input & Metadata page first:")
        show_issues(blockers)
        return
    if len(view.levels) < 2 or view.factor is None:
        st.warning("The comparison factor needs at least two levels in metadata.csv.")
        return
    unsafe = services.unsafe_levels(view.levels)
    if unsafe:
        st.error(f"Level(s) {', '.join(map(repr, unsafe))} cannot be used in a contrast ID; rename them in metadata.csv.")
        return
    for issue in report.warnings:
        if issue.code not in {"kegg_snapshot_not_pinned"}:
            st.warning(issue.message)

    existing = list(report.contrasts.contrasts) if report.contrasts else []
    if existing:
        st.subheader("Saved contrasts (contrasts.csv)")
        st.dataframe(pd.DataFrame([{"contrast_id": c.contrast_id, "numerator": c.numerator, "denominator": c.denominator, "direction": f"log2FC > 0 = higher in {c.numerator}"} for c in existing]), hide_index=True)

    st.subheader("Compose contrasts")
    key = f"contrasts::{project}"
    if key not in st.session_state:
        st.session_state[key] = [(c.numerator, c.denominator) for c in existing]
    chosen: list[tuple[str, str]] = st.session_state[key]
    stale = [pair for pair in chosen if pair[0] not in view.levels or pair[1] not in view.levels]
    if stale:
        # The metadata changed since these were composed; a contrast of a level
        # that no longer exists is never saved.
        chosen[:] = [pair for pair in chosen if pair not in stale]
        st.warning("Removed composed contrast(s) whose levels are no longer in metadata.csv: " + ", ".join(f"{num} vs {den}" for num, den in stale))
    saved_message = st.session_state.pop(f"{key}::saved", None)
    if saved_message:
        st.success(saved_message)
    cols = st.columns([2, 2, 1, 2])
    numerator = cols[0].selectbox("Numerator (test group)", view.levels, index=len(view.levels) - 1, key=f"{key}::num")
    denominator = cols[1].selectbox("Denominator (reference group)", view.levels, index=0, key=f"{key}::den")
    add_reason = ("Numerator and denominator must be different levels." if numerator == denominator
                  else "This contrast is already in the list." if (numerator, denominator) in chosen else None)
    if cols[2].button("Add", key=f"{key}::add", disabled=add_reason is not None, help=add_reason):
        chosen.append((numerator, denominator))
        st.rerun()
    if add_reason:
        st.caption(f"*Add* is disabled: {add_reason}")
    if cols[3].button("Propose all pairs (rnaseq design rule)", key=f"{key}::propose"):
        st.session_state[key] = [(item.numerator, item.denominator) for item in build_proposal(report, condition_order=view.levels).contrasts]
        st.rerun()
    if not chosen:
        st.info("No contrast selected yet.")
    for index, (num, den) in enumerate(chosen):
        definition = services.contrast_definition(view.factor, num, den)
        box = st.container(border=True)
        row = box.columns([5, 1])
        row[0].markdown(
            f"**{definition.contrast_id}** — numerator `{num}`, denominator `{den}`, formula `{view.formula}`  \n"
            f"Positive log2FC = higher in **{num}** (than in {den}); negative = higher in {den}."
        )
        if row[1].button("Remove", key=f"{key}::remove::{index}"):
            chosen.pop(index)
            st.rerun()

    selection = design_enrichment_selection(config, None)
    gene_id_type = None
    if design_needs_gene_id_type(config, selection):
        gene_id_type = st.selectbox("Gene identifier type in the count matrix (needed for enrichment)", [value for value, _ in wizard.GENE_ID_TYPE_CHOICES], key=f"{key}::geneid")
    if config.project.preset is Preset.L1:
        st.caption("Saving contrasts moves this project to L2, exactly as `rnaseq design` does.")
    replace = not existing or st.checkbox("Replace the saved contrasts", key=f"{key}::replace")
    save_reason = ("Add at least one contrast first." if not chosen
                   else "Tick *Replace the saved contrasts* to overwrite contrasts.csv." if not replace else None)
    if save_reason:
        st.caption(f"*Save contrasts* is disabled: {save_reason}")
    if st.button("Save contrasts", type="primary", key=f"{key}::save", disabled=save_reason is not None, help=save_reason):
        try:
            no_df = residual_df_blocker(report)
            if no_df is not None:
                raise DesignSetupError(no_df + " Add replicates or simplify the design before configuring L2.")
            definitions = [services.contrast_definition(view.factor, num, den) for num, den in chosen]
            selected, annotation = resolve_design_enrichment(config, None, gene_id_type)
            after = apply_design(report, services.design_proposal(report, view.levels, definitions), enrichment=selected, annotation=annotation)
        except DesignRejected as exc:
            st.error(str(exc))
            show_issues(exc.issues)
        except DesignSetupError as exc:
            st.error(str(exc))
        else:
            st.session_state[f"{key}::saved"] = f"Wrote {len(definitions)} contrast(s), froze the condition order ({', '.join(view.levels)}), scope L2; validation {after.state}. Planning artifacts refreshed."
            st.session_state.pop(f"{key}::replace", None)
            st.rerun()  # redraw the saved contrasts and the validation state from the files just written

    pin = services.kegg_pin_status(validated(project))
    if pin.applies:
        st.subheader("KEGG snapshot")
        if pin.pinned:
            st.success(f"Pinned to `{pin.pinned}` (annotation.enrichment.kegg.snapshot).")
        else:
            if pin.warning:
                st.warning(pin.warning)
            if pin.latest and st.button(f"Pin to the current snapshot {pin.latest}", key=f"{key}::pin"):
                try:
                    pin_kegg_snapshot(report.loaded.config_path, pin.latest)
                except (DesignSetupError, ValueError) as exc:
                    st.error(str(exc))
                else:
                    st.rerun()


# --------------------------------------------------------------------- E. Analysis runner

@st.cache_data(show_spinner="Checking the downstream runtime ...")
def _runtime_prefix(revision: str | None) -> tuple[str, bool, str | None]:
    from rnaseq.downstream_runtime import planned_runtime_prefix

    try:
        prefix, ready = planned_runtime_prefix()
    except Exception as exc:  # reported, never fatal for the page
        return "", False, str(exc)
    return str(prefix), ready, None


def _runtime_section(active: jobs.Launch | None) -> None:
    from rnaseq.downstream_runtime import runtime_source_revision

    if active is not None:
        st.caption("Runtime check skipped while a run is active.")
        return
    prefix, ready, error = _runtime_prefix(runtime_source_revision())
    if error:
        st.error(f"Downstream runtime check failed: {error}")
    elif ready:
        st.success(f"Downstream runtime ready: `{prefix}`")
    else:
        st.warning(f"The first run creates the downstream runtime (`{prefix}`); this takes several minutes before analysis starts.")


SCOPE_OUTPUTS = {
    Preset.QC: "FASTQ technical QC (fastp, MultiQC). No statistics.",
    Preset.L1: "L1: QC, filtering, DESeq2 normalization, blind VST, PCA, sample correlation and distribution plots. No differential expression.",
    Preset.L2: "L2: everything in L1, plus DESeq2 differential expression per contrast (raw LFC for DEG calls, shrunken LFC for display) and the selected enrichment analyses.",
}


def _summary(report, command: list[str]) -> None:
    config = report.config
    contrasts = report.contrasts.contrasts if report.contrasts else ()
    enrichment = tuple(config.analysis.enrichment) if config.analysis else ()
    rows: dict[str, object] = {
        "project": config.project.id,
        "scope": services.scope_label(report),
        "species": config.organism.species.value,
        "input": config.input.type.value + (f" — {len(report.counts.sample_ids)} samples, {report.counts.gene_count} genes" if report.counts else "") + (f" — {len(report.fastq.sample_ids)} samples, {config.input.layout.value}" if report.fastq else ""),
        "quantification": workspace.quantification_label(config),
        "strandedness": config.upstream.strandedness if config.input.type is InputType.FASTQ else None,
        "FASTQ preprocessing": (f"{config.input.preprocessing.value} — " + dict(wizard.PREPROCESSING_CHOICES).get(config.input.preprocessing.value, ""))
                               if config.input.type is InputType.FASTQ and config.input.preprocessing is not None else None,
        "reference": workspace.reference_label(config),
    }
    if config.project.preset is not Preset.QC:
        rows["design"] = f"{config.design.type.value}, {config.design.formula}"
    if config.project.preset is Preset.L2:
        rows["DEG thresholds"] = f"padj < {config.thresholds.padj}, |log2FC| ≥ {config.thresholds.abs_log2fc} (raw LFC)"
        for item in contrasts:
            rows[f"contrast {item.contrast_id}"] = f"log2FC > 0 = higher in {item.numerator} than {item.denominator}"
        rows["enrichment"] = describe_enrichment(enrichment)
        if config.annotation is not None and enrichment:
            rows["KEGG snapshot"] = config.annotation.enrichment.kegg.snapshot or "not pinned (latest local snapshot)"
    rows["resources"] = f"CPUs {config.execution.max_cpus}, memory GiB {config.execution.max_memory_gb}"
    facts_table(rows)
    st.caption(SCOPE_OUTPUTS[config.project.preset])
    if config.project.preset is Preset.L1:
        st.caption("To add differential expression (L2), define contrasts on the Contrasts page; saving contrasts moves the project to L2. "
                   "An L1 run can later be followed by an L2 run that reuses its upstream results.")
    with st.expander("Command"):
        st.code(shlex.join(command), language="bash", wrap_lines=True)


def _validation_overview(report) -> None:
    """The project validator's findings, grouped by workflow step, in plain words."""

    if report.is_valid:
        text = "The project is valid and can be run."
        if report.warnings:
            text += f" {len(report.warnings)} warning(s) below do not block the run."
        st.success(text)
    else:
        st.error(f"The project cannot run yet: {len(report.errors)} problem(s) to fix.")
    grouped: dict[str, list] = {}
    for issue in report.issues:
        grouped.setdefault(workspace.issue_area(issue), []).append(issue)
    for area in workspace.STEP_NAMES:
        issues = grouped.get(area)
        if not issues:
            continue
        target = STEP_PAGES.get(area)
        link = f" — fix on {PAGES[target].title}" if target in PAGES and target != "run" else ""
        st.markdown(f"**{area}**{link}")
        show_issues(issues)


STATUS_REFRESH_SECONDS = 10
LIVE_RECORDED_STATES = {"CREATED", "PLANNED", "RUNNING"}


def live_run_in_project(project: Path, active: jobs.Launch | None) -> bool:
    """Whether a run of this project is starting or running (cheap: run_state.json and process liveness only)."""

    if active is not None and active.project == str(project.resolve()):
        return True
    cache = _Cache()
    for run_dir, state in discover_runs(project):
        if str(state.get("status")) in LIVE_RECORDED_STATES:
            shown, _note = jobs.run_display_state(run_dir, state, cache)
            if workspace.execution_status(shown) == "RUNNING":
                return True
    return False


def _run_status_section(config: GuiConfig, project: Path, active: jobs.Launch | None) -> None:
    """Runs of the project; while one is active this section alone re-reads the run files every few seconds.

    It is a Streamlit fragment: an automatic refresh reruns only this
    function, so the project is not re-validated. When the run reaches a
    final state the whole page is redrawn once and the timer stops.
    """

    st.subheader("Runs of this project")
    live = live_run_in_project(project, active)
    st.fragment(_run_status_body, run_every=STATUS_REFRESH_SECONDS if live else None)(config, project, live)


def _run_status_body(config: GuiConfig, project: Path, live: bool) -> None:
    active = jobs.active_launch()
    if live != live_run_in_project(project, active):
        st.rerun()  # started or finished: redraw the page (Run button, timer) once
    cols = st.columns([1, 5])
    cols[0].button("Refresh", key="run_refresh")
    if live:
        cols[1].caption(f"A run of this project is active: this section updates automatically every {STATUS_REFRESH_SECONDS} s "
                        f"(last read {datetime.now().strftime('%H:%M:%S')}; the project is not re-validated). It stops when the run ends.")
    else:
        cols[1].caption(f"Read from the run files at {datetime.now().strftime('%H:%M:%S')}; press Refresh to update.")
    launches = jobs.list_launches(project)
    for launch in launches[:1]:
        run_dir = jobs.launch_run_dir(launch)
        alive, note = jobs.launch_alive(launch)
        if run_dir is None:
            if alive:
                st.info(f"Submission {launch.job_id} (case {launch.case_id}) is starting: `rnaseq run` is checking the project before it creates the run directory.")
            elif active is None or active.job_id != launch.job_id:
                tail = jobs.tail(launch.log_path, 40)
                if "Case run:" not in tail:
                    st.error(f"Submission {launch.job_id} (case {launch.case_id}) ended without creating a run ({note}). Launcher log `{launch.log}`:")
                    st.code(tail, language=None)
    cache = _Cache()
    runs = discover_runs(project)
    if not runs:
        st.write("No recorded runs.")
        return
    rows = workspace.project_runs(project, cache)
    st.dataframe(run_history_table(rows), hide_index=True, width="stretch")
    st.caption("Execution = whether the run finished (RUNNING, SUCCESS, FAILED, INCOMPLETE). Scientific QC = the QC verdicts the run recorded (PASS, WARNING, INCONCLUSIVE, NOT_EVALUATED); a QC warning does not change the execution status.")
    labels = [row.label for row in rows]
    chosen = remembered_choice("Run", labels, f"run_selected::{project}")
    row = rows[labels.index(chosen)]
    run_dir, state = runs[labels.index(chosen)]
    status = build_status(run_dir, state, cache)
    if row.execution == "INCOMPLETE":
        st.error(f"INCOMPLETE — the run did not finish; its files were not changed. {row.note or ''}")
    elif row.execution == "RUNNING" and state.get("status") == "CREATED":
        st.info("RUNNING — preparing: freezing inputs and creating or verifying the downstream runtime. With a new nf-rna build this first step takes several minutes; the analysis stages follow.")
    elif row.execution == "RUNNING":
        st.info("RUNNING — the stages below follow the run.")
    elif row.execution == "FAILED":
        st.error(f"FAILED: {state.get('error') or 'see the logs'}")
    elif row.execution == "SUCCESS":
        st.success(f"SUCCESS. Results: `{state.get('delivery') or run_dir / 'delivery'}`")
    if row.execution == "RUNNING":
        stop_instructions(run_dir, state)
    show_run_qc(row.qc)
    if row.has_results:
        cols = st.columns(2)
        if cols[0].button("Open results", key=f"run_open_results::{project}"):
            open_run(project, row, "results")
        if cols[1].button("Open report & downloads", key=f"run_open_reports::{project}"):
            open_run(project, row, "reports")
    with st.expander("Stages and tasks (rnaseq status)", expanded=row.execution == "RUNNING"):
        st.code(render_status(status, ascii_only=True), language=None)

    st.markdown("**Logs**")
    logs = sorted((run_dir / "logs").glob("*.log")) if (run_dir / "logs").is_dir() else []
    launch = next((item for item in launches if item.run_dir == str(run_dir)), None)
    options = {f"run: {path.name}": path for path in logs}
    if launch is not None:
        options[f"GUI launcher: {launch.log_path.name}"] = launch.log_path
    if not options:
        st.write("No logs yet.")
        return
    names = list(options)
    preferred = next((name for name in names if name == "run: rnaseq.log"), names[0])
    cols = st.columns([3, 1, 1])
    log_choice = cols[0].selectbox("Log", names, index=names.index(preferred), key=f"run_log::{project}::{chosen}")
    lines = cols[1].number_input("Last lines", min_value=10, max_value=5000, value=config.log_tail_lines, step=50, key="run_log_lines")
    alerts = cols[2].toggle("Only warnings and errors", key=f"run_log_alerts::{project}")
    st.caption(f"`{options[log_choice]}`")
    text = jobs.tail(options[log_choice], int(lines))
    if alerts:
        text = workspace.alert_lines(text) or "(no warning or error line in these lines)"
    st.code(text or "(empty)", language=None)


def project_fingerprint(report) -> tuple:
    """Modification times and sizes of the files a validation report depends on."""

    loaded = report.loaded
    if loaded is None:
        return ()
    paths = [loaded.config_path, loaded.metadata_path, loaded.contrasts_path, loaded.input_path]
    if loaded.input_path.is_dir():
        paths.extend(sorted(loaded.input_path.iterdir()))
    stamp = []
    for path in paths:
        try:
            info = path.stat()
        except OSError:
            stamp.append((str(path), None))
            continue
        stamp.append((str(path), info.st_mtime_ns, info.st_size))
    return tuple(stamp)


def run_blocked_reason(report, case_id: str, active: jobs.Launch | None) -> str | None:
    """Why *Run* cannot start now, in words (None when it can)."""

    from rnaseq.service import validate_case_id

    if active is not None:
        return f"a GUI run is active (case {active.case_id}); only one GUI run can be active at a time."
    if not report.is_valid:
        return f"the project has {len(report.errors)} validation problem(s) to fix (section 1)."
    if not case_id:
        return "enter a case ID (a name for this analysis, for example the date or the experiment)."
    try:
        validate_case_id(case_id)
    except (ExecutionPreflightError, ValueError) as exc:
        return f"the case ID is not usable: {exc}"
    return None


def stop_instructions(run_dir: Path, state: dict) -> None:
    """How to stop a running run by hand. There is no stop button (see docs/dev/gui_hardening_v170.md)."""

    record = results.read_json(run_dir / "logs" / "rnaseq.process.json") or {}
    pid = record.get("pid")
    with st.expander("How to stop this run"):
        st.markdown(
            "The GUI has no stop button: it cannot yet guarantee that every process of a run is stopped "
            "cleanly in every phase. To stop the run, use a terminal on the server."
        )
        if not isinstance(pid, int):
            st.caption("This run recorded no process identity (`logs/rnaseq.process.json`); find its `rnaseq run` process with `ps`.")
            return
        st.code(f"ps -o pid,lstart,args -p {pid}    # check that this is the rnaseq run of {run_dir.parent.name}\nkill -TERM {pid}", language="bash")
        st.markdown(
            "- `kill -TERM` (not `kill -9`) lets `rnaseq` stop Nextflow, which stops its tasks; the run is then recorded as "
            "**INTERRUPTED** and can be retried with `rnaseq retry`. `kill -9` leaves Nextflow and its tasks running.\n"
            "- Outside a Nextflow stage (preparing, or between stages) `rnaseq` stops at once without recording it; the GUI then shows the run as **INCOMPLETE**.\n"
            "- Do not stop a run that is **preparing** (recorded CREATED) while it creates the downstream runtime: the Conda "
            "environment would be left half-built, and later runs refuse it until it is removed by hand."
        )
        if state.get("status") == "CREATED":
            st.warning("This run is preparing (CREATED). Prefer to let it reach the analysis stages before stopping it.")


def runner_page() -> None:
    from rnaseq.cli import render_validation_report
    from rnaseq.planner import generate_plan

    config = gui_config()
    st.title("Validate & Run")
    project = require_project()
    if project is None:
        return
    report = validated(project)
    if report.config is None:
        show_issues(report.errors)
        return
    active = jobs.active_launch()
    st.caption(f"`{project}`")
    workflow_bar(project, report, "run")

    st.subheader("1. Validate")
    _validation_overview(report)
    fingerprint = project_fingerprint(report)
    st.caption("The findings above are re-read from the project files on every visit (read-only). The button below also writes the "
               "planning files under `planning/` (as `rnaseq plan`); it is optional: *Run* validates again before it starts.")
    if st.button("Validate and refresh the plan", key="run_validate"):
        text = render_validation_report(report)
        if report.is_valid:
            generate_plan(report)
            text += "\n\nPlanning artifacts refreshed (as `rnaseq plan`)."
        st.session_state[f"validation::{project}"] = (fingerprint, text)
    stored = st.session_state.get(f"validation::{project}")
    if stored is not None:
        if stored[0] != fingerprint:
            # Written for an earlier version of the project; never shown as current.
            st.session_state.pop(f"validation::{project}", None)
            st.caption("The full validation report was removed because project.yaml, the metadata, the contrasts or the input changed since it was made. The findings above are current.")
        else:
            with st.expander("Full validation report (rnaseq validate)", expanded=True):
                st.code(stored[1], language=None)

    st.subheader("2. Settings")
    preset = report.config.project.preset
    case_id = remembered_text("Case ID", f"case::{project}", placeholder=f"e.g. {preset.value}-{datetime.now():%Y%m%d}",
                              help="A name for this analysis: letters, digits, '.', '_' or '-'. Each press of Run creates a new immutable run "
                                   "under runs/<case ID>/<run ID> (the run ID is the start time), so a case ID can be reused for re-runs; "
                                   "nothing earlier is overwritten.")
    existing_cases = sorted({run_dir.parent.name for run_dir, _state in discover_runs(project)})
    if existing_cases:
        st.caption("Case IDs used before in this project: " + ", ".join(f"`{item}`" for item in existing_cases[:12]))
    reuse = None
    options = workspace.reuse_options(report)
    candidates = [item.reference for item in options if item.accepted]
    if len(candidates) == 1:
        reuse = candidates[0] if st.checkbox(f"Reuse the upstream results of {candidates[0]}", value=True, key=f"run_reuse_one::{project}") else None
    elif len(candidates) > 1:
        picked = st.selectbox("Reuse upstream results of", ["(do not reuse)", *candidates], key=f"run_reuse_many::{project}")
        reuse = None if picked == "(do not reuse)" else picked
    refused = [item for item in options if not item.accepted]
    if candidates:
        st.caption("Reuse skips FASTQ processing: the earlier run's quantification is used when its frozen upstream inputs (FASTQ, reference, settings) are identical.")
    if refused:
        with st.expander(f"Earlier runs that cannot be reused ({len(refused)})"):
            for item in refused:
                st.markdown(f"- `{item.reference}`: {item.reason}")
    elif report.config.input.type is InputType.RAW_COUNTS:
        st.caption("Raw-count projects have no upstream stage, so there is nothing to reuse.")
    command = jobs.rnaseq_command(project, case_id or "CASE-ID", reuse)
    _summary(report, command)
    if preset is not Preset.QC:
        _runtime_section(active)

    st.subheader("3. Run")
    label = {Preset.QC: "Run QC", Preset.L1: "Run L1", Preset.L2: "Run L2"}[preset]
    reason = run_blocked_reason(report, case_id, active)
    if active is not None:
        st.warning(f"A GUI run is active (case {active.case_id} of `{active.project}`); a new run can start when it ends.")
    elif reason:
        st.caption(f"*{label}* is disabled: {reason}")
    if st.button(label, type="primary", key="run_start", disabled=reason is not None, help=reason):
        try:
            with st.spinner("Checking the run preflight ..."):
                services.run_preflight(project, case_id, reuse)
            launch = jobs.launch_run(project, case_id, config, reuse_upstream=reuse)
        except jobs.ActiveRunError as exc:
            st.warning(str(exc))
        except (ExecutionPreflightError, jobs.LaunchError, ValueError) as exc:
            st.error(f"Not started: {exc}")
        except OSError as exc:
            st.error(f"System error: {exc}")
        else:
            st.success(f"Started case {case_id} ({launch.launcher} launcher). The run continues if this page or the GUI is closed.")
            active = launch
    _run_status_section(config, project, active)


# --------------------------------------------------------------------- F. Results explorer

@st.cache_data(max_entries=32, show_spinner=False)
def _table(path: str, mtime: float) -> results.Table:
    return results.read_table(Path(path))


def table(path: Path) -> results.Table:
    return _table(str(path), path.stat().st_mtime)


def frame(data: results.Table) -> pd.DataFrame:
    """A table for display: numeric columns as numbers (so sorting by a column header is numeric), the rest as text.

    Missing values (``NA``) become empty cells and sort last. Values are parsed
    from the file's text for display only; downloads are the files themselves.
    """

    table = pd.DataFrame(list(data.rows), columns=list(data.columns), dtype=str)
    for column, kind in results.numeric_columns(data).items():
        values = table[column].map(lambda text: None if (text or "").strip() in results.NA_TOKENS else float(text))
        table[column] = values.astype("Int64") if kind == "int" else values.astype("float64")
    return table


def show_frame(data: pd.DataFrame, **kwargs):
    return st.dataframe(data, column_config=number_format(data), **kwargs)


def number_format(data: pd.DataFrame) -> dict[str, object]:
    """Column settings that show floating-point statistics without rounding small p-values to 0."""

    return {str(column): st.column_config.NumberColumn(str(column), format="%.6g")
            for column, dtype in data.dtypes.items() if pd.api.types.is_float_dtype(dtype)}


def download(path: Path, label: str, key: str, config: GuiConfig) -> None:
    """One file as a download button, or its server path and scp command when it is too large."""

    if results.too_large(path, config.download_limit_mb):
        st.caption(f"{label}: larger than {config.download_limit_mb} MB; copy it from the server: `{path}`")
        st.code(results.scp_command(path), language="bash")
        return
    st.download_button(label, data=lambda: path.read_bytes(), file_name=path.name, key=key, on_click="ignore")


def tiff_download(option: figure_export.TiffOption, label: str, key: str, config: GuiConfig) -> None:
    """A TIFF download (the pipeline's file, or one encoded in memory on click), or a disabled button saying why not."""

    if not option.available:
        st.button(label, key=key, disabled=True, help=f"Unavailable: {option.reason}")
        return
    if option.method == figure_export.PIPELINE_FILE:
        if results.too_large(option.source, config.download_limit_mb):
            download(option.source, label, key, config)
            return
    st.download_button(label, data=lambda: figure_export.tiff_bytes(option), file_name=option.file_name, mime="image/tiff",
                       key=key, on_click="ignore", help=f"{option.dpi} dpi TIFF: {option.method}")


def png_download(option: figure_export.PngOption, label: str, key: str, config: GuiConfig) -> None:
    """A PNG at the submission dpi (the pipeline's file, or rendered from the PDF on click), or a disabled button saying why not."""

    if not option.available:
        st.button(label, key=key, disabled=True, help=f"Unavailable: {option.reason}")
        return
    if option.source.suffix == ".png":
        download(option.source, label, key, config)
        return
    st.download_button(label, data=lambda: figure_export.png_bytes(option), file_name=option.file_name, mime="image/png",
                       key=key, on_click="ignore", help=f"{option.dpi} dpi PNG: {option.method}")


def figure_downloads_menu(item: results.Figure, key: str, config: GuiConfig) -> None:
    """Every download of one figure in a single menu: vector PDF, PNG and TIFF at 300 dpi, TIFF at 600 dpi, and the PNG preview."""

    with st.popover("Download", icon=":material/download:"):
        st.caption(f"{item.label}. Nothing is redrawn from data: the pipeline's own files, or its vector PDF rendered at the stated dpi. "
                   "TIFF is RGB with LZW compression; a raster is never enlarged.")
        if item.pdf is not None:
            download(item.pdf, "PDF (vector)", f"{key}::pdf", config)
        else:
            st.button("PDF (vector)", key=f"{key}::pdf", disabled=True, help="Unavailable: the pipeline wrote no PDF for this figure")
        png_download(figure_export.png_option(item.stem, pdf=item.pdf), f"PNG {figure_export.PNG_DPI} dpi", f"{key}::png300", config)
        options = figure_export.tiff_options(item.stem, item.pdf)
        for dpi in figure_export.TIFF_DPIS:
            tiff_download(options[dpi], f"TIFF {dpi} dpi", f"{key}::tiff{dpi}", config)
        download(item.png, "PNG 150 dpi (report preview)", f"{key}::png", config)


def figure_caption(item: results.Figure) -> str:
    size = f" · {item.size_mm[0]:.0f} × {item.size_mm[1]:.0f} mm as drawn" if item.size_mm else ""
    return item.label + size


def show_figure(item: results.Figure, key: str, config: GuiConfig, *, width: str = "content") -> None:
    if item.png is None:
        st.warning(f"**{item.label}: unavailable** — {item.reason}")
        return
    st.image(str(item.png), caption=figure_caption(item), width=width)
    figure_downloads_menu(item, key, config)


# Figures without a recorded size are laid out as this wide (mm).
DEFAULT_FIGURE_MM = 90.0
# A row of figures is at most this wide (mm) before the next figure starts a new row.
FIGURE_ROW_MM = 220.0


def figure_rows(items: list[results.Figure]) -> list[list[results.Figure]]:
    """Figures packed into rows of at most FIGURE_ROW_MM physical width, in order."""

    rows: list[list[results.Figure]] = []
    width = 0.0
    for item in items:
        mm = item.size_mm[0] if item.size_mm else DEFAULT_FIGURE_MM
        if rows and width + mm <= FIGURE_ROW_MM:
            rows[-1].append(item)
            width += mm
        else:
            rows.append([item])
            width = mm
    return rows


def show_figure_grid(items: list[results.Figure], key: str, config: GuiConfig, scale_mm: float) -> None:
    """Figures side by side at one common on-screen scale.

    Column widths are proportional to each figure's physical width (and a
    spacer fills the row up to ``scale_mm``), so every figure on the page is
    shown at the same scale and its text at the same size, whatever the screen.
    """

    for row_index, row in enumerate(figure_rows(items)):
        widths = [item.size_mm[0] if item.size_mm else DEFAULT_FIGURE_MM for item in row]
        spacer = scale_mm - sum(widths)
        columns = st.columns([*widths, spacer] if spacer > 1 else widths, gap="medium")
        for index, (column, item) in enumerate(zip(columns, row)):
            with column:
                show_figure(item, f"{key}::{row_index}::{index}", config, width="stretch")


def facts_table(facts: dict[str, object]) -> None:
    rows = [{"item": name, "value": str(value)} for name, value in facts.items() if value not in (None, "")]
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")


def choose_run(project: Path, page: str) -> results.ResultRun | None:
    """The run shown on Results or Reports; the selection is kept per project (key ``<page>_run::<project>``)."""

    runs = results.result_runs(project, _Cache())
    if not runs:
        st.info("This project has no run with results yet; start one on the Run page.")
        return None
    labels = [run.label for run in runs]
    by_label = {run.label: run for run in runs}
    chosen = remembered_choice("Run", labels, f"{page}_run::{project}",
                               format_func=lambda label: f"{label} — {workspace.execution_status(by_label[label].state)}, {by_label[label].level or '?'}")
    run = by_label[chosen]
    execution = workspace.execution_status(run.state)
    if execution != "SUCCESS":
        st.warning(f"This run is {execution}; only the outputs it wrote are shown, and missing ones are marked unavailable.")
    return run


L1_GROUP_HELP = {
    results.L1_STRUCTURE: "PCA and Pearson correlation of the blind VST matrix (all retained genes): how samples group, and whether any sample stands apart.",
    results.L1_DISTRIBUTION: "Per-sample expression distributions before normalization (log2 raw count + 1) and after blind VST; distribution-level QC only.",
}


def _qc_section(run: results.ResultRun, config: GuiConfig) -> None:
    view = results.qc_view(run)
    summary = view.summary
    if summary is None:
        st.warning("**L1 QC: unavailable** — l1/backend_summary.json was not produced.")
    else:
        facts_table({
            "samples": len(summary.get("samples") or []), "genes input": summary.get("genes_input"),
            "all-zero genes removed": summary.get("genes_removed_all_zero"), "low-total genes removed": summary.get("genes_removed_low_total"),
            "genes retained": summary.get("genes_retained"), "filter": (summary.get("filter") or {}).get("rule"),
            "normalization": (summary.get("normalization") or {}).get("method"), "VST": (summary.get("vst") or {}).get("method"),
            "flags": ", ".join(map(str, summary.get("flags") or [])) or "none",
        })
    _strandedness_section(view)
    if view.multiqc is not None:
        st.markdown("**Upstream QC (MultiQC)**")
        download(view.multiqc, f"MultiQC report ({view.multiqc.name})", "qc::multiqc", config)
    groups: dict[str, list[results.Figure]] = {}
    for item in view.figures:
        groups.setdefault(item.group or "Figures", []).append(item)
    # One scale for every L1 figure: the widest row of any group sets it.
    scale = max([sum(item.size_mm[0] if item.size_mm else DEFAULT_FIGURE_MM for item in row) for items in groups.values() for row in figure_rows(items)] + [180.0])
    for name, items in groups.items():
        st.markdown(f"**{name}**")
        st.caption(L1_GROUP_HELP.get(name, ""))
        show_figure_grid(items, f"qc::{name}", config, scale)
    st.caption("All figures are shown at one common scale (their physical sizes are in the captions); open *Download* under a figure for "
               "the vector PDF, PNG 300 dpi or TIFF 300/600 dpi.")
    for name, path in view.tables.items():
        with st.expander(name):
            show_frame(frame(table(path)), hide_index=True, width="stretch")
            download(path, f"Download {name}", f"qc::table::{name}", config)


def _strandedness_section(view: results.QcView) -> None:
    """Strandedness consistency of FASTQ runs (report only; counts use the declared strandedness)."""

    strand = view.strandedness
    if not strand:
        return
    status = strand.get("status")
    samples = strand.get("samples") or {}
    if status == "WARNING":
        st.warning("**Strandedness check: WARNING** — the reads contradict the declared strandedness "
                   f"(`{strand.get('declared')}`) for at least one sample. The counts were made with the declared strandedness and "
                   "were not changed; check the library preparation protocol.")
    elif status == "ADVISORY":
        st.info("**Strandedness check: note** — declared unstranded, but the reads look stranded. Unstranded counting is valid; check the protocol.")
    elif status == "INCONCLUSIVE":
        st.info("**Strandedness check: inconclusive** for at least one sample (too few strand-assignable fragments or no clear direction).")
    elif status == "CONSISTENT":
        st.success(f"**Strandedness check:** consistent with the declared `{strand.get('declared')}` for every sample.")
    else:
        st.caption(f"Strandedness check: not evaluated — {strand.get('reason') or 'no strand evidence recorded'}.")
    if samples:
        st.dataframe(pd.DataFrame([
            {"sample": name, "assigned -s 1": item.get("assigned_forward_s1"), "assigned -s 2": item.get("assigned_reverse_s2"),
             "-s 1 / -s 2": item.get("forward_over_reverse"), "observed": item.get("observed"), "status": item.get("status"), "detail": item.get("message")}
            for name, item in samples.items() if isinstance(item, dict)
        ]), hide_index=True, width="stretch")


def _deg_section(run: results.ResultRun, contrast: str, config: GuiConfig) -> None:
    view = results.deg_view(run, contrast)
    cols = st.columns(5)
    for col, (label, name) in zip(cols, (("significant", "significant"), ("up", "upregulated"), ("down", "downregulated"))):
        value = view.counts[name]
        col.metric(f"{label} genes", "unavailable" if value is None else value)
    cols[3].metric("DEG threshold: padj <", str(view.padj_below))
    cols[4].metric("DEG threshold: |log2FC| ≥", str(view.abs_log2fc_at_least))
    st.caption(f"DEG rule (contrast_metadata.json): {view.rule}")
    shrinkage = view.shrinkage
    explanation = (
        "DEGs are called on the **raw** `log2FoldChange` (with padj); `log2FoldChange_shrunk` is the shrunken "
        "effect size, used only for display and ranking."
    )
    if not shrinkage:
        st.info(explanation + " Shrinkage: not recorded for this run (L2 results from before nf-rna v1.5).")
    elif shrinkage.get("status") == "success":
        st.info(explanation + f" Shrinkage method: **{shrinkage.get('method')}** (coef `{shrinkage.get('coef')}`), status **{shrinkage.get('status')}**.")
    else:
        st.error(explanation + f" Shrinkage **{shrinkage.get('status')}** (method {shrinkage.get('method')}): {shrinkage.get('error') or 'no error recorded'}. `log2FoldChange_shrunk` is not a valid shrunken estimate for this contrast.")
    warnings = view.r_warnings
    if warnings:
        with st.expander(f"R warnings recorded during the fit and this contrast ({sum(int(item.get('count') or 1) for item in warnings)}) — the analysis completed"):
            st.dataframe(pd.DataFrame([
                {"step": item.get("step"), "message": item.get("message"), "count": item.get("count"), "call": item.get("call")} for item in warnings
            ]), hide_index=True, width="stretch")
    figure_cols = st.columns(3)
    for index, item in enumerate(view.figures):
        with figure_cols[index]:
            show_figure(item, f"deg::{contrast}::{index}", config)

    st.markdown("**DEG table** (`all_genes.tsv`; values exactly as written by DESeq2/R)")
    path = view.root / "all_genes.tsv"
    if not path.is_file():
        st.warning("**all_genes.tsv: unavailable** — not produced for this contrast.")
        return
    genes = table(path)
    key = f"deg::{run.run_dir}::{contrast}"
    cols = st.columns([3, 2, 2, 2, 2, 1])
    query = cols[0].text_input("Search symbol or gene_id", key=f"{key}::query")
    only = cols[1].toggle("Only significant", key=f"{key}::only")
    padj_max = cols[2].number_input("padj ≤", min_value=0.0, max_value=1.0, value=None, format="%g", key=f"{key}::padj")
    lfc_min = cols[3].number_input("|log2FC| ≥", min_value=0.0, value=None, format="%g", key=f"{key}::lfc")
    sort_by = cols[4].selectbox("Sort by", [column for column in results.DEG_SORT_COLUMNS if column in genes.columns], key=f"{key}::sort")
    descending = cols[5].toggle("Desc.", key=f"{key}::desc")
    selected = results.filter_deg(genes, results.significant_ids(view), query=query, only_significant=only, padj_max=padj_max, abs_lfc_min=lfc_min, sort_by=sort_by, descending=descending)
    st.caption(f"{len(selected.rows)} of {len(genes.rows)} genes. 'Only significant' = the genes in significant.tsv.")
    first = [column for column in ("gene_id", "symbol", "log2FoldChange", "log2FoldChange_shrunk", "padj", "pvalue", "baseMean") if column in selected.columns]
    show_frame(frame(selected), hide_index=True, width="stretch", key=f"{key}::table", column_order=[*first, *(column for column in selected.columns if column not in first)])
    cols = st.columns(4)
    for col, name in zip(cols, ("all_genes", "significant", "upregulated", "downregulated")):
        if (view.root / f"{name}.tsv").is_file():
            with col:
                download(view.root / f"{name}.tsv", f"{name}.tsv", f"{key}::dl::{name}", config)


def _enrichment_section(run: results.ResultRun, contrast: str, config: GuiConfig) -> None:
    modules = results.enrichment_modules(run)
    if not modules:
        st.info("No functional enrichment was selected for this run (analysis.enrichment: []).")
        return
    labels = [module.label for module in modules]
    module = modules[labels.index(st.radio("Analysis", labels, horizontal=True, key=f"enr::{run.run_dir}::module"))]
    if module.summary is None:
        st.warning(f"**{module.label}: unavailable** — {module.reason}")
        return
    if module.legacy:
        st.warning(f"**Legacy results.** This {module.label} output was produced before nf-rna v1.5 and records no v1.5 significance semantics. Its significant terms follow that version's rule; they are shown as written and not mixed with v1.5 definitions.")
    facts_table(module.facts)
    sets = results.enrichment_results(module, contrast)
    names = [item.label for item in sets]
    result = sets[names.index(st.selectbox("Result set", names, key=f"enr::{run.run_dir}::{module.key}::set"))]
    key = f"enr::{run.run_dir}::{contrast}::{module.key}::{result.label}"
    status_text = f"Status **{result.status}**" + (f" — {result.reason}" if result.reason else "")
    (st.success if result.status == "SUCCESS" else st.info)(status_text)
    facts_table(result.facts)
    dotplots = [item for item in result.figures if "running score" not in item.label]
    running = [item for item in result.figures if "running score" in item.label]
    cols = st.columns(max(1, len(dotplots)))
    for index, item in enumerate(dotplots):
        with cols[index]:
            show_figure(item, f"{key}::fig::{index}", config)
    if running:
        with st.expander(f"GSEA running-score plots ({sum(1 for item in running if item.png)})", expanded=True):
            cols = st.columns(2)
            for index, item in enumerate(running):
                with cols[index % 2]:
                    show_figure(item, f"{key}::rs::{index}", config)

    show_all = st.toggle("Show all evaluated terms (default: significant only)", key=f"{key}::all")
    terms = result.all_terms if show_all else result.significant
    if terms is None:
        st.warning(f"**{'All-terms' if show_all else 'Significant-terms'} table: unavailable** — not produced in `{result.directory}`.")
        return
    found = table(terms)
    st.caption(f"`{terms.name}`: {len(found.rows)} terms. Select a row to see its member genes.")
    event = show_frame(frame(found), hide_index=True, width="stretch", on_select="rerun", selection_mode="single-row", key=f"{key}::terms")
    download(terms, f"Download {terms.name}", f"{key}::dl::terms", config)
    if result.members is None:
        st.warning("**Member genes: unavailable** — no member table in these results.")
        return
    members = table(result.members)
    id_column = results.term_id_column(members)
    options = list(dict.fromkeys(members.column(id_column))) if id_column else []
    picked = None
    rows = event.selection.rows if event is not None else []
    if rows and "ID" in found.columns:
        picked = found.rows[rows[0]]["ID"]
    choose = f"{key}::pathway"
    if picked is not None and st.session_state.get(f"{choose}::last") != picked:
        st.session_state[f"{choose}::last"] = picked
        if picked in options:
            st.session_state[choose] = picked
    if picked is not None and picked not in options:
        st.info(f"{picked}: member genes are recorded for significant pathways only.")
    pathway = st.selectbox("Pathway", options, index=None, placeholder="Select a row above or choose", key=choose)
    st.markdown(f"**Member genes** (`{result.members.name}`" + (", legacy member table" if result.members_legacy else "") + ")")
    if pathway is not None:
        genes = results.pathway_members(result.members, pathway)
        st.caption(f"{pathway}: {len(genes.rows)} member genes.")
        show_frame(frame(genes), hide_index=True, width="stretch", height=min(38 + 35 * len(genes.rows), 1200), key=f"{key}::members")
    download(result.members, f"Download {result.members.name}", f"{key}::dl::members", config)


def results_page() -> None:
    config = gui_config()
    st.title("Results")
    project = require_project()
    if project is None:
        return
    workflow_bar(project, validated(project), "results")
    run = choose_run(project, "results")
    if run is None:
        return
    st.caption(f"`{run.run_dir}` · analysis level {run.level or 'not recorded'}. Figures are the pipeline's own PNG files; each has a *Download* menu (PDF, PNG 300 dpi, TIFF 300/600 dpi). Nothing is recomputed.")
    sections = ["QC & exploratory", *(["Differential expression", "Enrichment"] if run.is_l2 else []), "Provenance"]
    section = st.radio("Section", sections, horizontal=True, key=f"results_section::{project}::{run.level}")
    contrast = None
    if section in {"Differential expression", "Enrichment"}:
        contrasts = results.contrast_ids(run)
        if not contrasts:
            st.warning("**L2 results: unavailable** — no contrast output under downstream/l2/contrasts.")
            return
        contrast = st.selectbox("Contrast", contrasts, key=f"results_contrast::{run.run_dir}")
    if section == "QC & exploratory":
        with st.expander(f"Scientific QC summary: {workspace.run_qc(run.run_dir).headline()}", expanded=workspace.run_qc(run.run_dir).level == "WARNING"):
            show_run_qc(workspace.run_qc(run.run_dir))
        _qc_section(run, config)
    elif section == "Differential expression":
        _deg_section(run, contrast, config)
    elif section == "Enrichment":
        _enrichment_section(run, contrast, config)
    else:
        st.dataframe(pd.DataFrame(results.provenance(run), columns=["item", "value"]), hide_index=True, width="stretch")


# --------------------------------------------------------------------- G. Reports & downloads

def reports_page() -> None:
    import streamlit.components.v1 as components

    config = gui_config()
    st.title("Reports & downloads")
    project = require_project()
    if project is None:
        return
    run = choose_run(project, "reports")
    if run is None:
        return
    st.subheader("Report")
    report = results.report_path(run)
    if report is None:
        st.warning("**report.html: unavailable** — the run has not produced its report.")
    else:
        st.caption(f"`{report}` (self-contained; figures are embedded).")
        download(report, f"Download {report.name}", "rep::report", config)
        if st.toggle("Show the report here", key=f"rep::show::{run.run_dir}"):
            components.html(report.read_text(encoding="utf-8"), height=900, scrolling=True)

    st.subheader("Files")
    items = results.downloads(run)
    groups = list(dict.fromkeys(item.group for item in items if item.group != "Report"))
    for group in groups:
        members = [item for item in items if item.group == group]
        cols = st.columns([4, 1])
        labels = [item.label for item in members]
        chosen = cols[0].selectbox(f"{group} ({len(members)})", labels, key=f"rep::{run.run_dir}::{group}")
        with cols[1]:
            st.write("")
            item = members[labels.index(chosen)]
            download(item.path, "Download", f"rep::{group}::dl", config)

    figure_downloads(run, config)

    st.subheader("Delivery package")
    size = results.delivery_size(run)
    if size is None:
        st.warning("**Delivery package: unavailable** — the run has no delivery/delivery_manifest.yaml.")
        return
    st.caption(f"`{run.delivery}` · {size / 1024 / 1024:.1f} MB in the delivery manifest.")
    if size > config.download_limit_mb * 1024 * 1024:
        st.info(f"Larger than the download limit ({config.download_limit_mb} MB, `download_limit_mb` in gui.yaml); copy the folder from the server:")
        st.code(results.scp_command(run.delivery), language="bash")
        return
    archive = results.delivery_archive_path(run)
    if st.button("Prepare delivery package (.zip)", key="rep::zip"):
        try:
            with st.spinner("Verifying the delivery manifest and writing the zip ..."):
                archive = results.build_delivery_archive(run)
        except OSError as exc:
            st.error(str(exc))
    if archive.is_file():
        download(archive, f"Download {archive.name}", "rep::zip::dl", config)
        st.caption(f"Written to `{archive}` after checking every file against delivery_manifest.yaml.")


FIGURE_FORMATS = ("PDF", *(f"TIFF {dpi} dpi" for dpi in figure_export.TIFF_DPIS), "PNG (report preview)")


def figure_downloads(run: results.ResultRun, config: GuiConfig) -> None:
    """One figure at a time in PDF, TIFF (300 dpi default, 600 dpi) or PNG; nothing is redrawn from data."""

    st.subheader("Figures (PDF / TIFF)")
    stems = results.figure_stems(run)
    if not stems:
        st.warning("**Figures: unavailable** — the run has no figure PDFs.")
        return
    base = run.run_dir / "downstream"
    labels = [str(stem.relative_to(base)) for stem in stems]
    cols = st.columns([3, 2, 1])
    chosen = cols[0].selectbox(f"Figure ({len(stems)})", labels, key=f"rep::{run.run_dir}::figure")
    form = cols[1].radio("Format", FIGURE_FORMATS, index=1, horizontal=True, key=f"rep::{run.run_dir}::figfmt")
    stem = stems[labels.index(chosen)]
    pdf, png = stem.with_suffix(".pdf"), stem.with_suffix(".png")
    with cols[2]:
        st.write("")
        if form == "PDF":
            download(pdf, "Download", "rep::fig::dl", config)
            note = "the pipeline's Cairo PDF, unchanged"
        elif form.startswith("TIFF"):
            option = figure_export.tiff_option(stem, int(form.split()[1]), pdf)
            tiff_download(option, "Download", "rep::fig::dl", config)
            note = option.method if option.available else f"unavailable — {option.reason}"
        elif png.is_file():
            download(png, "Download", "rep::fig::dl", config)
            note = "the pipeline's 150 dpi report preview, unchanged"
        else:
            st.button("Download", key="rep::fig::dl", disabled=True)
            note = "unavailable — the pipeline wrote no PNG preview for this figure"
    st.caption(f"{form}: {note}. TIFF is RGB, LZW, with the dpi in its resolution tags; a raster is never enlarged to reach a higher dpi.")


PAGES: dict[str, st.Page] = {}


def main() -> None:
    st.set_page_config(page_title="nf-rna", layout="wide")
    PAGES.update({
        "dashboard": st.Page(dashboard_page, title="Projects", url_path="dashboard", default=True),
        "datasets": st.Page(datasets_page, title="Datasets", url_path="datasets"),
        "project": st.Page(project_page, title="Project overview", url_path="project"),
        "new": st.Page(new_project_page, title="New Project", url_path="new"),
        "inputs": st.Page(inputs_page, title="Input & Metadata", url_path="inputs"),
        "contrasts": st.Page(contrasts_page, title="Contrasts", url_path="contrasts"),
        "run": st.Page(runner_page, title="Validate & Run", url_path="run"),
        "results": st.Page(results_page, title="Results", url_path="results"),
        "reports": st.Page(reports_page, title="Reports & Downloads", url_path="reports"),
    })
    page = st.navigation({
        "Workspace": [PAGES["dashboard"], PAGES["datasets"], PAGES["project"]],
        "Set up": [PAGES["new"], PAGES["inputs"], PAGES["contrasts"]],
        "Analysis": [PAGES["run"], PAGES["results"], PAGES["reports"]],
    }, position="sidebar")
    sidebar(gui_config())
    page.run()
