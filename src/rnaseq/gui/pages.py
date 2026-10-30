"""Streamlit pages: Dashboard, New Project, Input & Metadata, Contrasts, Run, Results, Reports & Downloads.

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
from rnaseq.gui import jobs, results, services
from rnaseq.gui.config import GuiConfig, PathNotAllowed, ensure_data_layout, gui_state_dir, load_config
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


# --------------------------------------------------------------------- shared context

def gui_config() -> GuiConfig:
    config = load_config()
    if not st.session_state.get("_layout_checked"):
        st.session_state["_layout_messages"] = [*config.notes, *ensure_data_layout(config)]
        st.session_state["_layout_checked"] = True
    return config


def current_project() -> Path | None:
    value = st.session_state.get(PROJECT_KEY) or st.query_params.get(PROJECT_KEY)
    if not value:
        return None
    path = Path(value)
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
        st.info("Select a project in the sidebar, or create one on the New Project page.")
    return project


def sidebar(config: GuiConfig) -> None:
    projects = services.discover_projects(config)
    current = current_project()
    options = [str(path) for path in projects]
    if current is not None and str(current) not in options:
        options.insert(0, str(current))
    with st.sidebar:
        st.caption(f"nf-rna {PIPELINE_VERSION}")
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


def path_picker(config: GuiConfig, label: str, key: str, *, want: str, suffixes: tuple[str, ...] = (), start: Path | None = None) -> Path | None:
    """A text field plus a directory listing; only paths inside the allowed roots are returned."""

    roots = config.available_roots()
    if not roots:
        st.error("No allowed root is available; check the GUI settings and mounted disks.")
        return None
    text_key, cwd_key = f"{key}_text", f"{key}_cwd"
    if cwd_key not in st.session_state:
        st.session_state[cwd_key] = str(start if start is not None and start.is_dir() else roots[0])
    typed = st.text_input(label, key=text_key, placeholder=f"A {'directory' if want == 'dir' else 'file'} inside {', '.join(map(str, roots))}")
    with st.expander("Browse", expanded=not typed):
        cwd = Path(st.session_state[cwd_key])
        top = st.columns([1, 3])
        root_choice = top[0].selectbox("Root", [str(root) for root in roots], key=f"{key}_root")
        if top[0].button("Go to root", key=f"{key}_goroot"):
            st.session_state[cwd_key] = root_choice
            st.rerun()
        top[1].markdown(f"Folder: `{cwd}`")
        try:
            entries = sorted(cwd.iterdir(), key=lambda item: item.name.lower())
        except OSError as exc:
            st.error(f"Cannot list {cwd}: {exc}")
            entries = []
        directories = [entry for entry in entries if entry.is_dir() and not entry.name.startswith(".")]
        files = [entry for entry in entries if entry.is_file() and not entry.name.startswith(".") and (not suffixes or entry.name.endswith(suffixes))]
        row = st.columns([3, 1, 1])
        sub = row[0].selectbox("Subfolders", [entry.name for entry in directories], index=None, placeholder="(none selected)", key=f"{key}_sub")
        if row[1].button("Open folder", key=f"{key}_open", disabled=sub is None):
            st.session_state[cwd_key] = str(cwd / sub)
            st.rerun()
        if row[2].button("Up", key=f"{key}_up", disabled=any(cwd == root for root in roots)):
            try:
                st.session_state[cwd_key] = str(config.resolve_allowed(cwd.parent))
            except PathNotAllowed:
                pass
            st.rerun()
        if want == "dir":
            if st.button("Use this folder", key=f"{key}_usedir"):
                st.session_state[text_key] = str(cwd)
                st.rerun()
        else:
            row = st.columns([3, 1])
            chosen = row[0].selectbox("Files", [entry.name for entry in files], index=None, placeholder="(none selected)", key=f"{key}_file")
            if row[1].button("Use this file", key=f"{key}_usefile", disabled=chosen is None):
                st.session_state[text_key] = str(cwd / chosen)
                st.rerun()
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


def dashboard_page() -> None:
    config = gui_config()
    st.title("nf-rna")
    st.caption(f"Version {PIPELINE_VERSION} · settings: `{config.source}`")
    for message in st.session_state.get("_layout_messages", []):
        st.info(message)
    active = jobs.active_launch()
    if active is not None:
        st.warning(f"Active GUI run: case **{active.case_id}** of `{active.project}` (submitted {active.launched_at}). Open the Run page of that project to follow it.")

    st.subheader("Storage")
    st.dataframe(pd.DataFrame([{"allowed root": str(item.path), "status": "available" if item.available else "UNAVAILABLE", "detail": item.message} for item in config.root_statuses()]), hide_index=True, width="stretch")

    st.subheader("Projects")
    cache = _Cache()
    projects = services.discover_projects(config)
    rows = [services.project_row(path, cache) for path in projects]
    if rows:
        st.dataframe(pd.DataFrame([{"project": row.project_id, "scope": row.preset, "input": row.input_type, "latest run": row.latest_run, "state": row.latest_state, "path": str(row.path), "problem": row.error or ""} for row in rows]), hide_index=True, width="stretch")
        chosen = st.selectbox("Open project", [str(row.path) for row in rows], index=None, format_func=lambda value: Path(value).name, key="dash_open")
        if chosen and st.button("Open", key="dash_open_btn"):
            select_project(Path(chosen))
            if "inputs" in PAGES:
                st.switch_page(PAGES["inputs"])
            st.rerun()
    else:
        st.write("No projects found under the allowed roots.")

    st.subheader("Recent runs")
    recent = services.recent_runs(projects, cache=cache)
    if recent:
        st.dataframe(pd.DataFrame(recent), hide_index=True, width="stretch")
    else:
        st.write("No recorded runs.")

    left, right = st.columns(2)
    with left:
        st.subheader("Runtime (rnaseq doctor)")
        if st.button("Re-run checks", key="doctor_refresh"):
            _doctor.clear()
        checks = _doctor()
        st.dataframe(pd.DataFrame(checks, columns=["check", "verdict", "detail"]), hide_index=True, width="stretch")
    with right:
        st.subheader("Registered references")
        references, error = services.registered_references()
        if error:
            st.error(error)
        if references:
            st.dataframe(pd.DataFrame([{"species": ref.species, "provider": ref.provider, "release": ref.release, "assembly": ref.assembly_identity, "Salmon": "yes" if ref.salmon_index else "no", "HISAT2": "yes" if ref.hisat2_index else "no", "root": str(ref.root)} for ref in references]), hide_index=True, width="stretch")
        elif not error:
            st.write("No registered references.")


# --------------------------------------------------------------------- B. New Project

def new_project_page() -> None:
    from rnaseq.execution import detect_local_resource_capacity, validate_local_execution_budget
    from rnaseq.project import create_project
    from rnaseq.references import compatible_registered_references
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
        method = cols[0].selectbox("Quantification", [value for value, _ in wizard.METHOD_CHOICES], format_func=lambda value: dict(wizard.METHOD_CHOICES)[value], key="np_method")
        layout = cols[1].selectbox("Sequencing layout", [value for value, _ in wizard.LAYOUT_CHOICES], format_func=lambda value: f"{value} ({dict(wizard.LAYOUT_CHOICES)[value]})", key="np_layout")
        strand_options = [value for value, _ in wizard.strandedness_choices(method)]
        if method == "hisat2_featurecounts":
            strandedness = cols[2].selectbox("Strandedness", strand_options, index=None, placeholder="Choose explicitly", key="np_strand_hisat2")
            st.info(STRANDEDNESS_HELP)
        else:
            strandedness = cols[2].selectbox("Strandedness", strand_options, index=strand_options.index(wizard.default_strandedness(method)), format_func=lambda value: f"{value} ({dict(wizard.strandedness_choices(method))[value]})", key="np_strand_salmon")
        try:
            registered = compatible_registered_references(wizard.species_value(species).value, method)
        except Exception as exc:  # a broken registry must not block the other reference routes
            registered = []
            st.warning(f"Reference registry could not be read: {exc}")
        labels = {f"{ref.provider} {ref.release} / {ref.assembly_identity} ({ref.root})": ref for ref in registered}
        other = ["iGenomes (Salmon only)"] if method == "salmon" else []
        choice = st.selectbox("Reference", [*labels, *other, "Other (advanced)"], key="np_reference")
        with st.expander("Advanced settings"):
            preprocessing = st.selectbox("FASTQ preprocessing", [value for value, _ in wizard.PREPROCESSING_CHOICES], format_func=lambda value: f"{value} — {dict(wizard.PREPROCESSING_CHOICES)[value]}", key="np_pre")
            if choice == "Other (advanced)":
                source = st.selectbox("Reference source", ["local", "custom"], key="np_ref_source")
                if source == "local":
                    root = st.text_input("Managed reference root", key="np_ref_root")
                    manifest = st.text_input("Reference manifest", value="reference_manifest.yaml", key="np_ref_manifest")
                    reference_inputs = {"source": "local", "local_root": Path(root) if root else None, "local_manifest": manifest}
                else:
                    reference_inputs = {"source": "custom", "fasta": st.text_input("Genome FASTA", key="np_fasta") or None, "gtf": st.text_input("Annotation GTF", key="np_gtf") or None, "hisat2_index": st.text_input("Prepared HISAT2 index (optional)", key="np_hisat2") or None}
            else:
                reference_inputs = None
        try:
            if choice in labels:
                root, manifest = wizard.managed_reference_settings(labels[choice])
                reference = wizard.reference_options(source="local", species=wizard.species_value(species), method=method, local_root=root, local_manifest=manifest, fasta=None, gtf=None, transcript_fasta=None, salmon_index=None, hisat2_index=None)
            elif choice in other:
                reference = wizard.reference_options(source="igenomes", species=wizard.species_value(species), method=method, local_root=None, local_manifest=None, fasta=None, gtf=None, transcript_fasta=None, salmon_index=None, hisat2_index=None)
            elif reference_inputs is not None:
                to_path = lambda value: Path(value) if value else None  # noqa: E731
                reference = wizard.reference_options(
                    source=str(reference_inputs["source"]), species=wizard.species_value(species), method=method,
                    local_root=reference_inputs.get("local_root"), local_manifest=reference_inputs.get("local_manifest"),
                    fasta=to_path(reference_inputs.get("fasta")), gtf=to_path(reference_inputs.get("gtf")), transcript_fasta=None,
                    salmon_index=None, hisat2_index=to_path(reference_inputs.get("hisat2_index")),
                )
        except ProjectCreationError as exc:
            st.error(str(exc))

    preset = st.radio("Analysis scope", [value for value, _ in wizard.preset_choices(input_type)], format_func=lambda value: f"{value} — {dict(wizard.preset_choices(input_type))[value]}", key=f"np_preset_{input_type}", index=1 if input_type == "fastq" else 0)
    design_type, pair_id, batch, enrichment, gene_id_type = "two_group", None, False, None, None
    if preset != "qc":
        design_type = st.selectbox("Experimental design", [value for value, _ in wizard.DESIGN_CHOICES], format_func=lambda value: f"{value} — {dict(wizard.DESIGN_CHOICES)[value]}", key="np_design")
        if design_type == "paired_two_group":
            pair_id = st.text_input("Pairing metadata column", value=wizard.DEFAULT_PAIR_ID, key="np_pair", help="The metadata column naming the biological pairing unit (patient, donor, cell line ...).")
        batch = st.checkbox("There is a known technical batch (adds a 'batch' column; the model adjusts for it)", key="np_batch")
        if preset == "L2":
            enrichment = st.checkbox(f"Run functional enrichment: {describe_enrichment(DEFAULT_L2_ENRICHMENT)}", value=True, key="np_enrichment")
        if preset == "L1" or enrichment:
            gene_id_type = st.selectbox("Gene identifier type in the count matrix", [value for value, _ in wizard.GENE_ID_TYPE_CHOICES], format_func=lambda value: f"{value} — {dict(wizard.GENE_ID_TYPE_CHOICES)[value]}", key="np_geneid")

    suggestion = wizard.resource_suggestion()
    with st.expander("Local execution resources"):
        st.caption(f"Detected {suggestion.logical_cpus or '?'} CPUs / {suggestion.memory_gib or '?'} GiB; `auto` gives {suggestion.auto_cpus} CPUs / {suggestion.auto_memory_gib} GiB on this machine (resolved again on each run).")
        cpus = st.number_input("CPU limit (0 = auto)", min_value=0, value=0, step=1, key="np_cpus")
        memory = st.number_input("Memory limit in GiB (0 = auto)", min_value=0, value=0, step=1, key="np_memory")

    if st.button("Create project", type="primary", key="np_create"):
        try:
            parent = config.resolve_allowed(parent_text)
            if not name:
                raise ProjectCreationError("Enter a project name.")
            if input_type == "fastq" and method == "hisat2_featurecounts" and not strandedness:
                raise ProjectCreationError("Choose the strandedness explicitly for HISAT2 + featureCounts.")
            if input_type == "fastq" and reference is None:
                raise ProjectCreationError("Choose a reference.")
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
            st.page_link(PAGES["inputs"], label="Next: Input & Metadata →")


# --------------------------------------------------------------------- C. Input & Metadata

def _input_section(config: GuiConfig, project: Path, report) -> None:
    from rnaseq.project import fastq_samplesheet_text, import_count_matrix, import_fastq_samplesheet

    assert report.config is not None and report.loaded is not None
    st.subheader("1. Input data")
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
    else:
        st.warning(f"No FASTQ files in `{report.loaded.input_path}` yet.")
    mode = st.radio("Add FASTQ files", ["Detect in a folder", "Import an existing samplesheet"], horizontal=True, key="fq_mode")
    if mode == "Detect in a folder":
        folder = path_picker(config, "FASTQ folder", "fastq_dir", want="dir", start=config.project_parent.parent / "fastq")
        if folder is None:
            return
        summary, errors = services.discover_fastq(folder, layout)
        show_issues(errors)
        if summary is None or not summary.records:
            return
        st.dataframe(pd.DataFrame([{"sample_id": record.sample_id, "lane": record.lane, "R1": record.fastq_1.name, "R2": record.fastq_2.name if record.fastq_2 else ""} for record in summary.records]), hide_index=True, width="stretch")
        sheet = fastq_samplesheet_text(summary.records, report.config.upstream.strandedness or "auto")
        st.code(sheet, language=None)
        st.download_button("Download samplesheet.csv", sheet, file_name="samplesheet.csv", key="fq_download")
        if not errors and st.button(f"Copy these {len(summary.records)} FASTQ assignment(s) into the project", key="fq_import"):
            # The import keeps its own copy as planning/imported_fastq_samplesheet.csv.
            scratch = gui_state_dir() / "samplesheets" / f"{project.name}-{datetime.now():%Y%m%d-%H%M%S}.csv"
            scratch.parent.mkdir(parents=True, exist_ok=True)
            scratch.write_text(sheet, encoding="utf-8")
            try:
                _import_samplesheet(project, scratch, import_fastq_samplesheet)
            finally:
                scratch.unlink(missing_ok=True)
    else:
        sheet_path = path_picker(config, "Samplesheet CSV (sample,fastq_1,fastq_2,strandedness)", "fq_sheet", want="file", suffixes=(".csv",))
        if sheet_path is not None and st.button("Import samplesheet (copies the FASTQ files)", key="fq_sheet_import"):
            _import_samplesheet(project, sheet_path, import_fastq_samplesheet)


def _import_samplesheet(project: Path, sheet: Path, importer) -> None:
    try:
        with st.spinner("Copying FASTQ files into the project ..."):
            records = importer(project, sheet)
    except (ProjectCreationError, ProjectConfigError) as exc:
        st.error(str(exc))
        return
    except OSError as exc:
        st.error(f"System error while copying: {exc}")
        return
    st.success(f"Imported {len(records)} samplesheet row(s).")
    st.rerun()


def _metadata_state_key(project: Path) -> str:
    return f"metadata::{project}"


def _set_editor_table(project: Path, frame: pd.DataFrame) -> None:
    st.session_state[_metadata_state_key(project)] = frame
    st.session_state[f"{_metadata_state_key(project)}::version"] = st.session_state.get(f"{_metadata_state_key(project)}::version", 0) + 1


def _metadata_section(project: Path, report) -> None:
    assert report.config is not None
    st.subheader("2. Sample metadata")
    if report.config.project.preset is Preset.QC:
        st.info("QC-only projects make no statistical design claim; metadata is not used.")
        return
    columns = services.metadata_columns(report)
    design_vars = services.design_variables(report)
    st.caption(
        f"Design formula `{report.config.design.formula}`; required columns: sample_id, {', '.join(design_vars)}. "
        "Conditions are entered by you; they are never derived from file names."
    )
    key = _metadata_state_key(project)
    if key not in st.session_state:
        rows = services.merge_sample_ids(services.read_metadata(report), columns, services.input_sample_ids(report))
        _set_editor_table(project, pd.DataFrame(rows, columns=columns))
    frame: pd.DataFrame = st.session_state[key]

    tools = st.tabs(["Paste table (TSV / Excel)", "Import CSV", "Assign values to samples", "Reload"])
    with tools[0]:
        pasted = st.text_area("Paste rows copied from Excel (tab-separated) or CSV, header row first", key=f"{key}::paste", height=150)
        if st.button("Apply pasted table", key=f"{key}::paste_apply"):
            header, rows = services.parse_table_text(pasted)
            if "sample_id" not in header:
                st.error("The pasted table needs a header row with a sample_id column.")
            else:
                merged_columns = list(dict.fromkeys([*columns, *header]))
                _set_editor_table(project, pd.DataFrame(services.merge_sample_ids(rows, merged_columns, services.input_sample_ids(report)), columns=merged_columns))
                st.rerun()
    with tools[1]:
        uploaded = st.file_uploader("Metadata CSV/TSV", type=["csv", "tsv", "txt"], key=f"{key}::upload")
        if uploaded is not None and st.button("Load into the editor", key=f"{key}::upload_apply"):
            header, rows = services.parse_table_text(uploaded.getvalue().decode("utf-8-sig"))
            if "sample_id" not in header:
                st.error("The file needs a sample_id column.")
            else:
                merged_columns = list(dict.fromkeys([*columns, *header]))
                _set_editor_table(project, pd.DataFrame(services.merge_sample_ids(rows, merged_columns, []), columns=merged_columns))
                st.rerun()
    with tools[2]:
        # A fresh widget key after each Apply clears the previous selection.
        bulk = f"{key}::bulk::{st.session_state.get(f'{key}::version', 0)}"
        cols = st.columns([3, 2, 2, 1])
        samples = cols[0].multiselect("Samples", [str(value) for value in frame["sample_id"].tolist() if str(value)], key=f"{bulk}::samples")
        column = cols[1].selectbox("Column", [column for column in frame.columns if column != "sample_id"], key=f"{bulk}::column")
        value = cols[2].text_input("Value", key=f"{bulk}::value")
        if cols[3].button("Apply", key=f"{bulk}::apply", disabled=not samples or column is None):
            updated = frame.copy()
            updated.loc[updated["sample_id"].astype(str).isin(samples), column] = value
            _set_editor_table(project, updated)
            st.rerun()
    with tools[3]:
        if st.button("Discard unsaved edits and reload metadata.csv", key=f"{key}::reload"):
            st.session_state.pop(key, None)
            st.rerun()

    edited = st.data_editor(
        frame, num_rows="dynamic", hide_index=True, width="stretch",
        key=f"{key}::editor::{st.session_state.get(f'{key}::version', 0)}",
        column_config={column: st.column_config.TextColumn(column, required=column in ("sample_id", *design_vars)) for column in frame.columns},
    )
    if st.button("Save metadata", type="primary", key=f"{key}::save"):
        rows = edited.fillna("").astype(str).to_dict("records")
        result = services.save_metadata(project, list(edited.columns), rows)
        st.session_state[f"{key}::result"] = result
        if result.saved:
            st.session_state[key] = edited
    result = st.session_state.get(f"{key}::result")
    if result is not None:
        if result.saved:
            st.success(f"Saved {report.loaded.metadata_path}.")
        else:
            st.error("metadata.csv was not changed; fix these first:")
        by_column: dict[str, list] = {}
        for issue in [*result.errors, *result.incomplete, *result.agreement]:
            by_column.setdefault(services.issue_column(issue, list(edited.columns)), []).append(issue)
        for column, issues in by_column.items():
            st.markdown(f"**{column}**")
            show_issues(issues)
        if result.saved and not (result.incomplete or result.agreement):
            st.caption("All design values are filled. Next: Contrasts (L2) or Run.")


def inputs_page() -> None:
    from rnaseq.validators import validate_project

    config = gui_config()
    st.title("Input & metadata")
    project = require_project()
    if project is None:
        return
    report = validate_project(project)
    if report.config is None:
        show_issues(report.errors)
        return
    st.caption(f"`{project}` · {services.scope_label(report)} · input {report.config.input.type.value}")
    _input_section(config, project, report)
    _metadata_section(project, report)


# --------------------------------------------------------------------- D. Contrast designer

def contrasts_page() -> None:
    from rnaseq.design_setup import (
        DesignRejected, DesignSetupError, apply_design, build_proposal, design_blockers, design_enrichment_selection,
        design_needs_gene_id_type, pin_kegg_snapshot, residual_df_blocker, resolve_design_enrichment,
    )
    from rnaseq.validators import validate_project

    gui_config()
    st.title("Contrasts")
    project = require_project()
    if project is None:
        return
    report = validate_project(project)
    config = report.config
    if config is None:
        show_issues(report.errors)
        return
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
    cols = st.columns([2, 2, 1, 2])
    numerator = cols[0].selectbox("Numerator (test group)", view.levels, index=len(view.levels) - 1, key=f"{key}::num")
    denominator = cols[1].selectbox("Denominator (reference group)", view.levels, index=0, key=f"{key}::den")
    if cols[2].button("Add", key=f"{key}::add", disabled=numerator == denominator):
        if (numerator, denominator) not in chosen:
            chosen.append((numerator, denominator))
        st.rerun()
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
    if st.button("Save contrasts", type="primary", key=f"{key}::save", disabled=not chosen or not replace):
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
            st.success(f"Wrote {len(definitions)} contrast(s), froze the condition order ({', '.join(view.levels)}), scope L2; validation {after.state}. Planning artifacts refreshed.")

    pin = services.kegg_pin_status(validate_project(project))
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


def _summary(report, command: list[str]) -> None:
    config = report.config
    contrasts = report.contrasts.contrasts if report.contrasts else ()
    enrichment = tuple(config.analysis.enrichment) if config.analysis else ()
    lines = [
        f"- Project: **{config.project.id}** · scope {services.scope_label(report)}",
        f"- Input: {config.input.type.value}" + (f" ({len(report.counts.sample_ids)} samples, {report.counts.gene_count} genes)" if report.counts else "") + (f" ({len(report.fastq.sample_ids)} samples)" if report.fastq else ""),
    ]
    if config.project.preset is not Preset.QC:
        lines.append(f"- Design: {config.design.type.value}, `{config.design.formula}`")
        lines.extend(f"- Contrast {c.contrast_id}: log2FC > 0 = higher in **{c.numerator}** than {c.denominator}" for c in contrasts)
    if config.project.preset is Preset.L2:
        lines.append(f"- Enrichment: {describe_enrichment(enrichment)}")
        if config.annotation is not None and enrichment:
            lines.append(f"- KEGG snapshot: {config.annotation.enrichment.kegg.snapshot or 'not pinned (latest local snapshot)'}")
    lines.append(f"- Resources: CPUs {config.execution.max_cpus}, memory GiB {config.execution.max_memory_gb}")
    st.markdown("\n".join(lines))
    st.markdown("Command:")
    st.code(shlex.join(command), language="bash", wrap_lines=True)


def _run_status_section(config: GuiConfig, project: Path, active: jobs.Launch | None) -> None:
    st.subheader("Runs")
    if st.button("Refresh", key="run_refresh"):
        st.rerun()
    st.caption(f"Read from the run files at {datetime.now().strftime('%H:%M:%S')}; refresh manually.")
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
    table = []
    for run_dir, state in runs:
        shown, note = jobs.run_display_state(run_dir, state, cache)
        table.append({"case": run_dir.parent.name, "run": run_dir.name, "state": shown, "started": state.get("started_at"), "completed": state.get("completed_at") or "", "note": note or ""})
    st.dataframe(pd.DataFrame(table), hide_index=True, width="stretch")
    labels = [f"{row['case']}/{row['run']}" for row in table]
    chosen = st.selectbox("Run", labels, index=0, key="run_selected")
    run_dir, state = runs[labels.index(chosen)]
    shown, note = jobs.run_display_state(run_dir, state, cache)
    status = build_status(run_dir, state, cache)
    if shown == "INCOMPLETE":
        st.error(f"INCOMPLETE — the run's process no longer exists; the run did not finish. The run files were not changed. {note}")
    elif shown == "RUNNING" and state.get("status") == "CREATED":
        st.info("RUNNING — preparing: freezing inputs and creating or verifying the downstream runtime. With a new nf-rna build this first step takes several minutes; the analysis stages follow.")
    elif shown == "FAILED":
        st.error(f"FAILED: {state.get('error') or 'see the logs'}")
    elif shown == "SUCCESS":
        st.success(f"SUCCESS. Results: `{state.get('delivery') or run_dir / 'delivery'}`")
        if "results" in PAGES:
            st.page_link(PAGES["results"], label="Open in Results →")
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
    cols = st.columns([3, 1])
    log_choice = cols[0].selectbox("Log", list(options), key="run_log")
    lines = cols[1].number_input("Last lines", min_value=10, max_value=5000, value=config.log_tail_lines, step=50, key="run_log_lines")
    st.caption(f"`{options[log_choice]}`")
    st.code(jobs.tail(options[log_choice], int(lines)) or "(empty)", language=None)


def runner_page() -> None:
    from rnaseq.cli import render_validation_report
    from rnaseq.planner import generate_plan
    from rnaseq.validators import validate_project

    config = gui_config()
    st.title("Run")
    project = require_project()
    if project is None:
        return
    report = validate_project(project)
    if report.config is None:
        show_issues(report.errors)
        return
    active = jobs.active_launch()

    st.subheader("Validate")
    if st.button("Validate", key="run_validate"):
        text = render_validation_report(report)
        if report.is_valid:
            generate_plan(report)
            text += "\n\nPlanning artifacts refreshed (as `rnaseq plan`)."
        st.session_state[f"validation::{project}"] = text
    if f"validation::{project}" in st.session_state:
        st.code(st.session_state[f"validation::{project}"], language=None)

    st.subheader("Settings")
    preset = report.config.project.preset
    case_id = st.text_input("Case ID", key=f"case::{project}", placeholder=f"e.g. {preset.value}-{datetime.now():%Y%m%d}", help="A new immutable run is created under runs/<case ID>/<run ID>.")
    reuse = None
    candidates = services.reuse_candidates(report)
    if len(candidates) == 1:
        reuse = candidates[0] if st.checkbox(f"Reuse the upstream results of {candidates[0]}", value=True, key="run_reuse_one") else None
    elif len(candidates) > 1:
        picked = st.selectbox("Reuse upstream results of", ["(do not reuse)", *candidates], key="run_reuse_many")
        reuse = None if picked == "(do not reuse)" else picked
    command = jobs.rnaseq_command(project, case_id or "CASE-ID", reuse)
    _summary(report, command)
    if preset is not Preset.QC:
        _runtime_section(active)

    label = {Preset.QC: "Run QC", Preset.L1: "Run L1", Preset.L2: "Run L2"}[preset]
    if active is not None:
        st.warning(f"A GUI run is active (case {active.case_id} of `{active.project}`); a new run can start when it ends.")
    if st.button(label, type="primary", key="run_start", disabled=active is not None or not case_id):
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
    return pd.DataFrame(list(data.rows), columns=list(data.columns), dtype=str)


def download(path: Path, label: str, key: str, config: GuiConfig) -> None:
    """One file as a download button, or its server path and scp command when it is too large."""

    if results.too_large(path, config.download_limit_mb):
        st.caption(f"{label}: larger than {config.download_limit_mb} MB; copy it from the server: `{path}`")
        st.code(results.scp_command(path), language="bash")
        return
    st.download_button(label, data=lambda: path.read_bytes(), file_name=path.name, key=key, on_click="ignore")


def show_figure(item: results.Figure, key: str, config: GuiConfig) -> None:
    if item.png is None:
        st.warning(f"**{item.label}: unavailable** — {item.reason}")
        return
    st.image(str(item.png), caption=item.label)
    cols = st.columns(2)
    with cols[0]:
        download(item.png, "PNG", f"{key}::png", config)
    if item.pdf is not None:
        with cols[1]:
            download(item.pdf, "PDF", f"{key}::pdf", config)


def facts_table(facts: dict[str, object]) -> None:
    rows = [{"item": name, "value": str(value)} for name, value in facts.items() if value not in (None, "")]
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")


def choose_run(project: Path, key: str) -> results.ResultRun | None:
    runs = results.result_runs(project, _Cache())
    if not runs:
        st.info("This project has no run with results yet; start one on the Run page.")
        return None
    labels = [f"{run.label} — {run.state}, {run.level or '?'}" for run in runs]
    chosen = st.selectbox("Run", labels, index=0, key=key)
    run = runs[labels.index(chosen)]
    if run.state != "SUCCESS":
        st.warning(f"This run is {run.state}; only the outputs it wrote are shown, and missing ones are marked unavailable.")
    return run


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
    if view.multiqc is not None:
        st.markdown("**Upstream QC (MultiQC)**")
        download(view.multiqc, f"MultiQC report ({view.multiqc.name})", "qc::multiqc", config)
    cols = st.columns(2)
    for index, item in enumerate(view.figures):
        with cols[index % 2]:
            show_figure(item, f"qc::{index}", config)
    for name, path in view.tables.items():
        with st.expander(name):
            st.dataframe(frame(table(path)), hide_index=True, width="stretch")
            download(path, f"Download {name}", f"qc::table::{name}", config)


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
    key = f"deg::{run.label}::{contrast}"
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
    st.dataframe(frame(selected), hide_index=True, width="stretch", key=f"{key}::table", column_order=[*first, *(column for column in selected.columns if column not in first)])
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
    module = modules[labels.index(st.radio("Analysis", labels, horizontal=True, key=f"enr::{run.label}::module"))]
    if module.summary is None:
        st.warning(f"**{module.label}: unavailable** — {module.reason}")
        return
    if module.legacy:
        st.warning(f"**Legacy results.** This {module.label} output was produced before nf-rna v1.5 and records no v1.5 significance semantics. Its significant terms follow that version's rule; they are shown as written and not mixed with v1.5 definitions.")
    facts_table(module.facts)
    sets = results.enrichment_results(module, contrast)
    names = [item.label for item in sets]
    result = sets[names.index(st.selectbox("Result set", names, key=f"enr::{run.label}::{module.key}::set"))]
    key = f"enr::{run.label}::{contrast}::{module.key}::{result.label}"
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
    event = st.dataframe(frame(found), hide_index=True, width="stretch", on_select="rerun", selection_mode="single-row", key=f"{key}::terms")
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
        st.dataframe(frame(genes), hide_index=True, width="stretch", height=min(38 + 35 * len(genes.rows), 1200), key=f"{key}::members")
    download(result.members, f"Download {result.members.name}", f"{key}::dl::members", config)


def results_page() -> None:
    config = gui_config()
    st.title("Results")
    project = require_project()
    if project is None:
        return
    run = choose_run(project, "results_run")
    if run is None:
        return
    st.caption(f"`{run.run_dir}` · analysis level {run.level or 'not recorded'}. Figures are the pipeline's own PNG files (PDF for download); nothing is recomputed.")
    sections = ["QC & exploratory", *(["Differential expression", "Enrichment"] if run.is_l2 else []), "Provenance"]
    section = st.radio("Section", sections, horizontal=True, key=f"results_section::{run.level}")
    contrast = None
    if section in {"Differential expression", "Enrichment"}:
        contrasts = results.contrast_ids(run)
        if not contrasts:
            st.warning("**L2 results: unavailable** — no contrast output under downstream/l2/contrasts.")
            return
        contrast = st.selectbox("Contrast", contrasts, key=f"results_contrast::{run.label}")
    if section == "QC & exploratory":
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
    run = choose_run(project, "reports_run")
    if run is None:
        return
    st.subheader("Report")
    report = results.report_path(run)
    if report is None:
        st.warning("**report.html: unavailable** — the run has not produced its report.")
    else:
        st.caption(f"`{report}` (self-contained; figures are embedded).")
        download(report, f"Download {report.name}", "rep::report", config)
        if st.toggle("Show the report here", key="rep::show"):
            components.html(report.read_text(encoding="utf-8"), height=900, scrolling=True)

    st.subheader("Files")
    items = results.downloads(run)
    groups = list(dict.fromkeys(item.group for item in items if item.group != "Report"))
    for group in groups:
        members = [item for item in items if item.group == group]
        cols = st.columns([4, 1])
        labels = [item.label for item in members]
        chosen = cols[0].selectbox(f"{group} ({len(members)})", labels, key=f"rep::{group}")
        with cols[1]:
            st.write("")
            item = members[labels.index(chosen)]
            download(item.path, "Download", f"rep::{group}::dl", config)

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


PAGES: dict[str, st.Page] = {}


def main() -> None:
    st.set_page_config(page_title="nf-rna", layout="wide")
    PAGES.update({
        "dashboard": st.Page(dashboard_page, title="Dashboard", url_path="dashboard", default=True),
        "new": st.Page(new_project_page, title="New Project", url_path="new"),
        "inputs": st.Page(inputs_page, title="Input & Metadata", url_path="inputs"),
        "contrasts": st.Page(contrasts_page, title="Contrasts", url_path="contrasts"),
        "run": st.Page(runner_page, title="Run", url_path="run"),
        "results": st.Page(results_page, title="Results", url_path="results"),
        "reports": st.Page(reports_page, title="Reports & Downloads", url_path="reports"),
    })
    page = st.navigation(list(PAGES.values()), position="sidebar")
    sidebar(gui_config())
    page.run()
