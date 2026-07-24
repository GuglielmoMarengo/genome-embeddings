"""Genome Embeddings scientific dashboard.

Install dependencies with ``python -m pip install -r requirements.txt`` and run
with ``python app.py``.  The mathematical workflow lives in ``src.dashboard``;
this module is only the NiceGUI presentation layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import plotly.graph_objects as go

try:
    from nicegui import events, ui
except ModuleNotFoundError:  # Allows core tests before optional UI install.
    events = None
    ui = None

from src.dashboard import (
    DashboardAnalysis,
    DashboardConfig,
    DatasetRecord,
    analyze_records,
    load_demo_records,
    parse_fasta_bytes,
)
from src.genome import GenomeMatrix


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_K_VALUES = [1, 2, 3, 4, 5]


@dataclass(slots=True)
class DashboardState:
    records: list[DatasetRecord] = field(default_factory=list)
    analysis: DashboardAnalysis | None = None


APP_CSS = """
:root {
  --ge-ink: #152033;
  --ge-muted: #667085;
  --ge-surface: #ffffff;
  --ge-background: #f4f7fb;
  --ge-primary: #2563eb;
  --ge-secondary: #0f766e;
  --ge-border: #dce3ee;
}
body { background: var(--ge-background); color: var(--ge-ink); }
.ge-shell { max-width: 1500px; margin: 0 auto; width: 100%; }
.ge-card {
  background: var(--ge-surface);
  border: 1px solid var(--ge-border);
  border-radius: 18px;
  box-shadow: 0 8px 30px rgba(30, 45, 75, 0.06);
}
.ge-kpi {
  min-height: 122px;
  background: linear-gradient(145deg, #ffffff, #f8fbff);
}
.ge-label { color: var(--ge-muted); font-size: 0.78rem; letter-spacing: .06em; text-transform: uppercase; }
.ge-value { color: var(--ge-ink); font-size: 1.7rem; font-weight: 700; }
.ge-section-title { font-size: 1.15rem; font-weight: 700; }
.ge-subtle { color: var(--ge-muted); }
.ge-brand-mark {
  width: 42px; height: 42px; border-radius: 13px;
  display: grid; place-items: center; color: white; font-weight: 800;
  background: linear-gradient(135deg, #2563eb, #0f766e);
  box-shadow: 0 8px 20px rgba(37, 99, 235, .25);
}
.q-tab--active { font-weight: 700; }
"""


def matrix_figure(matrix: GenomeMatrix) -> go.Figure:
    figure = go.Figure(
        data=go.Heatmap(
            z=matrix.values,
            x=matrix.labels,
            y=matrix.labels,
            colorscale="Viridis",
            colorbar={"title": metric_title(matrix.metric)},
            hovertemplate=(
                "%{y} vs %{x}<br>Value: %{z:.6f}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"{metric_title(matrix.metric)} · k={matrix.kmer_length}",
        margin={"l": 120, "r": 40, "t": 70, "b": 110},
        height=610,
        template="plotly_white",
        xaxis={"tickangle": -35},
        yaxis={"autorange": "reversed"},
    )
    return figure


def trajectory_figure(
    trajectories: dict[str, dict[int, float]],
    title: str,
    y_axis_title: str,
) -> go.Figure:
    figure = go.Figure()
    for name, trajectory in trajectories.items():
        figure.add_trace(
            go.Scatter(
                x=list(trajectory),
                y=list(trajectory.values()),
                mode="lines+markers",
                name=name,
                hovertemplate="k=%{x}<br>%{y:.6f}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        title=title,
        xaxis_title="k-mer length",
        yaxis_title=y_axis_title,
        template="plotly_white",
        height=460,
        legend={"orientation": "h", "y": 1.12},
        margin={"l": 70, "r": 30, "t": 90, "b": 60},
    )
    return figure


def ranking_stability_figure(analysis: DashboardAnalysis) -> go.Figure:
    stability = analysis.jensen_shannon_ranking_stability
    transitions = [f"{first}->{second}" for first, second in stability]
    values = [float(row["kendall_tau"]) for row in stability.values()]
    figure = go.Figure(
        go.Scatter(
            x=transitions,
            y=values,
            mode="lines+markers",
            name="Jensen–Shannon",
            hovertemplate="%{x}<br>Kendall τ=%{y:.3f}<extra></extra>",
        )
    )
    figure.add_hline(y=1.0, line_dash="dash", opacity=0.45)
    figure.add_hline(y=0.0, line_dash="dot", opacity=0.35)
    figure.update_layout(
        title="Jensen–Shannon ranking stability",
        xaxis_title="Scale transition",
        yaxis_title="Kendall rank correlation",
        yaxis={"range": [-1.05, 1.05]},
        template="plotly_white",
        height=420,
        margin={"l": 70, "r": 30, "t": 75, "b": 60},
    )
    return figure


def metric_title(metric: str) -> str:
    return {
        "euclidean": "Legacy Euclidean distance",
        "cosine": "Legacy cosine similarity",
        "euclidean_v2": "Descriptor V2 Euclidean distance",
        "cosine_v2": "Descriptor V2 cosine similarity",
        "embedding_v2_euclidean": "Multiscale embedding V2 distance",
        "embedding_v2_cosine": "Multiscale embedding V2 cosine similarity",
        "jensen_shannon": "Jensen–Shannon distance",
    }[metric]


def unique_label(records: list[DatasetRecord], requested: str) -> str:
    existing = {record.label for record in records}
    if requested not in existing:
        return requested
    index = 2
    while f"{requested} ({index})" in existing:
        index += 1
    return f"{requested} ({index})"


def create_dashboard() -> None:
    if ui is None:
        raise RuntimeError(
            "NiceGUI is not installed. Run: "
            "python -m pip install -r requirements.txt"
        )

    state = DashboardState(records=load_demo_records(PROJECT_ROOT))

    ui.add_head_html(f"<style>{APP_CSS}</style>")
    ui.colors(primary="#2563eb", secondary="#0f766e", accent="#7c3aed")

    with ui.header().classes("bg-white text-slate-900 border-b border-slate-200"):
        with ui.row().classes("ge-shell items-center justify-between px-5 py-3"):
            with ui.row().classes("items-center gap-3"):
                ui.html('<div class="ge-brand-mark">GE</div>')
                with ui.column().classes("gap-0"):
                    ui.label("Genome Embeddings").classes("text-xl font-bold")
                    ui.label("Turning genomes into mathematics.").classes(
                        "text-xs text-slate-500"
                    )
            dark_mode = ui.dark_mode()
            ui.button(
                icon="dark_mode",
                on_click=dark_mode.toggle,
            ).props("flat round").tooltip("Toggle dark mode")

    with ui.column().classes("ge-shell p-5 gap-5"):
        with ui.card().classes("ge-card w-full p-5"):
            with ui.row().classes("w-full items-start gap-5"):
                with ui.column().classes("grow gap-3"):
                    ui.label("Analysis workspace").classes("ge-section-title")
                    ui.label(
                        "Load single-record FASTA files or use the demonstration "
                        "dataset, choose k-mer scales, and compare the legacy, V2, "
                        "multiscale embedding, and Jensen–Shannon representations."
                    ).classes("ge-subtle max-w-4xl")
                    upload = ui.upload(
                        label="Add FASTA files",
                        multiple=True,
                        auto_upload=True,
                    ).props('accept=".fasta,.fa,.fna,.txt" flat bordered').classes(
                        "w-full"
                    )
                with ui.column().classes("w-80 gap-3"):
                    k_values_select = ui.select(
                        options=list(range(1, 9)),
                        value=DEFAULT_K_VALUES,
                        label="k-mer scales",
                        multiple=True,
                    ).props("outlined use-chips").classes("w-full")
                    selected_k_select = ui.select(
                        options=DEFAULT_K_VALUES,
                        value=3,
                        label="Single-scale k",
                    ).props("outlined").classes("w-full")
                    reference_select = ui.select(
                        options=[record.label for record in state.records],
                        value=state.records[0].label,
                        label="Reference sequence",
                    ).props("outlined").classes("w-full")
                    comparison_select = ui.select(
                        options=[record.label for record in state.records],
                        value=state.records[1].label,
                        label="Comparison sequence",
                    ).props("outlined").classes("w-full")
                    with ui.row().classes("w-full"):
                        analyze_button = ui.button("Run analysis", icon="science")
                        reset_button = ui.button("Reset demo", icon="restart_alt").props(
                            "outline"
                        )

        status = ui.label().classes("text-sm text-slate-500 px-1")

        tabs = ui.tabs().classes("w-full bg-white rounded-xl border border-slate-200")
        with tabs:
            overview_tab = ui.tab("Overview", icon="dashboard")
            dataset_tab = ui.tab("Dataset", icon="dataset")
            descriptors_tab = ui.tab("Descriptors V2", icon="functions")
            matrices_tab = ui.tab("Matrices", icon="grid_on")
            multiscale_tab = ui.tab("Multiscale", icon="timeline")
            exports_tab = ui.tab("Exports", icon="download")

        with ui.tab_panels(tabs, value=overview_tab).classes(
            "w-full bg-transparent p-0"
        ):
            with ui.tab_panel(overview_tab).classes("p-0"):
                overview_container = ui.column().classes("w-full gap-5")
            with ui.tab_panel(dataset_tab).classes("p-0"):
                dataset_container = ui.column().classes("w-full gap-5")
            with ui.tab_panel(descriptors_tab).classes("p-0"):
                descriptors_container = ui.column().classes("w-full gap-5")
            with ui.tab_panel(matrices_tab).classes("p-0"):
                matrices_container = ui.column().classes("w-full gap-5")
            with ui.tab_panel(multiscale_tab).classes("p-0"):
                multiscale_container = ui.column().classes("w-full gap-5")
            with ui.tab_panel(exports_tab).classes("p-0"):
                exports_container = ui.column().classes("w-full gap-5")

    def refresh_selects() -> None:
        labels = [record.label for record in state.records]
        for selector in (reference_select, comparison_select):
            selector.options = labels
            if selector.value not in labels:
                selector.value = labels[0] if labels else None
            selector.update()
        if (
            reference_select.value == comparison_select.value
            and len(labels) > 1
        ):
            comparison_select.value = labels[1]
            comparison_select.update()

    def render_analysis() -> None:
        analysis = state.analysis
        if analysis is None:
            return
        summary = analysis.summary()

        overview_container.clear()
        with overview_container:
            with ui.row().classes("w-full gap-4"):
                kpis = [
                    ("Sequences", summary["dataset_size"], "dataset"),
                    ("Selected k", summary["selected_k"], "tag"),
                    (
                        "Descriptor V2 distance",
                        f"{summary['descriptor_v2_euclidean_distance']:.5f}",
                        "functions",
                    ),
                    (
                        "Jensen–Shannon distance",
                        f"{summary['jensen_shannon_distance']:.5f}",
                        "scatter_plot",
                    ),
                ]
                for label, value, icon in kpis:
                    with ui.card().classes("ge-card ge-kpi grow p-4"):
                        ui.icon(icon).classes("text-blue-600 text-2xl")
                        ui.label(label).classes("ge-label")
                        ui.label(str(value)).classes("ge-value")
            with ui.card().classes("ge-card w-full p-5"):
                ui.label("Representation comparison").classes("ge-section-title")
                rows = [
                    {
                        "representation": key.replace("_", " "),
                        "value": f"{float(value):.6f}",
                    }
                    for key, value in summary.items()
                    if key.endswith(("distance", "similarity"))
                ]
                ui.table(
                    columns=[
                        {
                            "name": "representation",
                            "label": "Representation",
                            "field": "representation",
                            "align": "left",
                        },
                        {
                            "name": "value",
                            "label": "Pair value",
                            "field": "value",
                        },
                    ],
                    rows=rows,
                    row_key="representation",
                ).classes("w-full")

        dataset_container.clear()
        with dataset_container:
            with ui.card().classes("ge-card w-full p-5"):
                ui.label("Dataset inventory").classes("ge-section-title")
                dataset_rows = [
                    {
                        **row,
                        "gc_content": f"{float(row['gc_content']) * 100:.2f}%",
                    }
                    for row in analysis.dataset_rows
                ]
                ui.table(
                    columns=[
                        {"name": "label", "label": "Label", "field": "label", "align": "left"},
                        {"name": "source", "label": "Source", "field": "source", "align": "left"},
                        {"name": "length", "label": "Length", "field": "length"},
                        {"name": "gc_content", "label": "GC", "field": "gc_content"},
                        {"name": "header", "label": "FASTA header", "field": "header", "align": "left"},
                    ],
                    rows=dataset_rows,
                    row_key="label",
                    pagination=10,
                ).classes("w-full")

        descriptors_container.clear()
        with descriptors_container:
            with ui.card().classes("ge-card w-full p-5"):
                ui.label("Finite-sample and sparsity diagnostics").classes(
                    "ge-section-title"
                )
                rows = []
                for row in analysis.descriptor_v2_rows:
                    rows.append(
                        {
                            "label": row["label"],
                            "length": row["length"],
                            "conditional_entropy": f"{row['conditional_entropy']:.4f}",
                            "finite_sample_entropy": f"{row['finite_sample_entropy']:.4f}",
                            "effective_kmer_count": f"{row['effective_kmer_count']:.2f}",
                            "theoretical_coverage": f"{row['theoretical_coverage']:.4f}",
                            "observable_coverage": f"{row['observable_coverage']:.4f}",
                            "singleton_fraction": f"{row['singleton_fraction']:.4f}",
                            "repeated_window_fraction": f"{row['repeated_window_fraction']:.4f}",
                        }
                    )
                ui.table(
                    columns=[
                        {"name": "label", "label": "Sequence", "field": "label", "align": "left"},
                        {"name": "length", "label": "Length", "field": "length"},
                        {"name": "conditional_entropy", "label": "Conditional H", "field": "conditional_entropy"},
                        {"name": "finite_sample_entropy", "label": "Finite-sample H", "field": "finite_sample_entropy"},
                        {"name": "effective_kmer_count", "label": "Effective k-mers", "field": "effective_kmer_count"},
                        {"name": "theoretical_coverage", "label": "Theoretical coverage", "field": "theoretical_coverage"},
                        {"name": "observable_coverage", "label": "Observable coverage", "field": "observable_coverage"},
                        {"name": "singleton_fraction", "label": "Singleton fraction", "field": "singleton_fraction"},
                        {"name": "repeated_window_fraction", "label": "Repeated fraction", "field": "repeated_window_fraction"},
                    ],
                    rows=rows,
                    row_key="label",
                    pagination=10,
                ).classes("w-full")

        matrices_container.clear()
        with matrices_container:
            matrix_items = [
                analysis.legacy_euclidean_matrix,
                analysis.descriptor_v2_euclidean_matrix,
                analysis.embedding_v2_euclidean_matrix,
                analysis.jensen_shannon_matrix,
                analysis.legacy_cosine_matrix,
                analysis.descriptor_v2_cosine_matrix,
            ]
            for matrix in matrix_items:
                with ui.card().classes("ge-card w-full p-3"):
                    ui.plotly(matrix_figure(matrix)).classes("w-full")

        multiscale_container.clear()
        with multiscale_container:
            with ui.card().classes("ge-card w-full p-3"):
                ui.plotly(
                    trajectory_figure(
                        {
                            "Legacy Euclidean": (
                                analysis.legacy_euclidean_pair_trajectory
                            ),
                            "Descriptor V2 Euclidean": (
                                analysis.descriptor_v2_euclidean_pair_trajectory
                            ),
                            "Jensen–Shannon": (
                                analysis.jensen_shannon_pair_trajectory
                            ),
                        },
                        (
                            f"{analysis.config.reference_label} vs "
                            f"{analysis.config.comparison_label}"
                        ),
                        "Pair value",
                    )
                ).classes("w-full")
            with ui.card().classes("ge-card w-full p-3"):
                ui.plotly(ranking_stability_figure(analysis)).classes("w-full")
            with ui.card().classes("ge-card w-full p-5"):
                ui.label("Jensen–Shannon ranking trajectory").classes(
                    "ge-section-title"
                )
                rows = [
                    {
                        "k": k,
                        "ranking": " > ".join(label for label, _ in ranking),
                    }
                    for k, ranking in (
                        analysis.jensen_shannon_ranking_trajectory.items()
                    )
                ]
                ui.table(
                    columns=[
                        {"name": "k", "label": "k", "field": "k"},
                        {"name": "ranking", "label": "Ranking", "field": "ranking", "align": "left"},
                    ],
                    rows=rows,
                    row_key="k",
                ).classes("w-full")

        exports_container.clear()
        with exports_container:
            with ui.card().classes("ge-card w-full p-5"):
                ui.label("Export analysis").classes("ge-section-title")
                ui.label(
                    "Exports are generated in memory. Uploaded filenames are not "
                    "used as filesystem paths."
                ).classes("ge-subtle")
                with ui.row().classes("gap-3 mt-3"):
                    ui.button(
                        "Download JSON",
                        icon="data_object",
                        on_click=lambda: ui.download(
                            analysis.to_json().encode("utf-8"),
                            "genome_embeddings_analysis.json",
                        ),
                    )
                    ui.button(
                        "Download summary CSV",
                        icon="table_view",
                        on_click=lambda: ui.download(
                            analysis.summary_csv().encode("utf-8-sig"),
                            "genome_embeddings_summary.csv",
                        ),
                    ).props("outline")

    def run_analysis() -> None:
        try:
            raw_k_values = k_values_select.value or []
            k_values = tuple(sorted(int(value) for value in raw_k_values))
            if selected_k_select.value not in k_values:
                selected_k_select.value = k_values[0] if k_values else None
                selected_k_select.update()
            config = DashboardConfig(
                k_values=k_values,
                selected_k=int(selected_k_select.value),
                reference_label=str(reference_select.value),
                comparison_label=str(comparison_select.value),
            )
            analyze_button.props("loading")
            state.analysis = analyze_records(state.records, config)
            render_analysis()
            status.text = (
                f"Analysis complete · {len(state.records)} sequences · "
                f"k={list(k_values)}"
            )
            ui.notify("Analysis complete", type="positive")
        except Exception as error:  # UI boundary: show validated core error.
            status.text = f"Analysis failed: {error}"
            ui.notify(str(error), type="negative", multi_line=True)
        finally:
            analyze_button.props(remove="loading")

    async def handle_upload(event: events.UploadEventArguments) -> None:
        try:
            data = await event.file.read()
            provisional = parse_fasta_bytes(event.file.name, data)
            record = DatasetRecord(
                label=unique_label(state.records, provisional.label),
                genome=provisional.genome,
                source=provisional.source,
            )
            state.records.append(record)
            refresh_selects()
            status.text = f"Added {record.label} ({record.genome.length()} nt)."
            ui.notify(f"Added {record.label}", type="positive")
        except Exception as error:
            ui.notify(str(error), type="negative", multi_line=True)

    def reset_demo() -> None:
        state.records = load_demo_records(PROJECT_ROOT)
        refresh_selects()
        run_analysis()
        upload.reset()

    def update_selected_k_options() -> None:
        values = sorted(int(value) for value in (k_values_select.value or []))
        selected_k_select.options = values
        if selected_k_select.value not in values:
            selected_k_select.value = values[0] if values else None
        selected_k_select.update()

    upload.on_upload(handle_upload)
    k_values_select.on_value_change(lambda _: update_selected_k_options())
    analyze_button.on_click(run_analysis)
    reset_button.on_click(reset_demo)

    run_analysis()


def main() -> None:
    if ui is None:
        raise SystemExit(
            "NiceGUI is not installed. Run: "
            "python -m pip install -r requirements.txt"
        )
    create_dashboard()
    ui.run(
        title="Genome Embeddings",
        favicon="🧬",
        reload=False,
        language="en-US",
    )


if __name__ in {"__main__", "__mp_main__"}:
    main()
