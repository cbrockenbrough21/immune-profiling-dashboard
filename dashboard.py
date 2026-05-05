from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

# ---------------------------------------------------------------------------
# Path constants
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
OUTPUTS_DIR = BASE_DIR / "outputs"
RESPONDER_DIR = OUTPUTS_DIR / "responder_analysis"

FREQ_CSV = OUTPUTS_DIR / "sample_population_frequencies.csv"

RESPONDER_CSV_TMPL = RESPONDER_DIR / "responder_differential_analysis_time_{label}.csv"
RESPONDER_PNG_TMPL = RESPONDER_DIR / "responder_stratification_boxplot_time_{label}.png"

BASELINE_CSV = OUTPUTS_DIR / "baseline_cohort_summary.csv"
LONGITUDINAL_CSV = OUTPUTS_DIR / "longitudinal_trajectories.csv"
OVERVIEW_CSV = OUTPUTS_DIR / "study_overview.csv"
BREAKDOWN_CSV = OUTPUTS_DIR / "study_breakdown.csv"

# ---------------------------------------------------------------------------
# Streamlit page config
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Immune Trial Analysis", layout="wide")
st.title("Immune Trial Analysis")

# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


@st.cache_data(show_spinner=False)
def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8").strip()


def _safe_exists(path: Path) -> bool:
    return path.exists() and path.is_file()


def _load_responder_df(label: str) -> pd.DataFrame | None:
    path = Path(str(RESPONDER_CSV_TMPL).format(label=label))
    if not _safe_exists(path):
        return None
    df = _read_csv(path)
    for c in ["p_value", "q_value"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _responder_headline(df: pd.DataFrame, label: str) -> tuple[str, str]:
    """
    Returns (status_line, detail_line) for the selected timepoint.
    """
    sig = df.loc[df["q_value"] < 0.05, "population"].tolist() if "q_value" in df.columns else []
    if sig:
        status = "Statistically significant differences detected (FDR 0.05)."
        detail = "Significant populations: " + ", ".join(sig)
        return status, detail

    status = "No immune populations meet the FDR 0.05 significance threshold."
    detail = ""

    if "q_value" in df.columns and df["q_value"].notna().any():
        best_row = df.loc[df["q_value"].idxmin()]
        best_pop = str(best_row["population"])
        best_q = float(best_row["q_value"])
        if best_q < 0.10:
            if label == "all":
                detail = (
                    f"Closest signal in pooled analysis: {best_pop} (q = {best_q:.3f}). "
                    "This is a trend, not statistically significant under the predefined threshold."
                )
            else:
                detail = (
                    f"Closest signal at day {label}: {best_pop} (q = {best_q:.3f}). "
                    "This is a trend, not statistically significant under the predefined threshold."
                )

    return status, detail


# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tab_overview, tab_frequencies, tab_responder, tab_baseline, tab_longitudinal = st.tabs(
    ["Overview", "Population Frequencies", "Responder Analysis", "Baseline Cohort", "Longitudinal Trajectories"]
)

# ---------------------------------------------------------------------------
# Overview
# ---------------------------------------------------------------------------
with tab_overview:
    st.markdown(
        "A multi-cohort clinical trial profiling peripheral blood immune cell populations "
        "across melanoma and NSCLC patients treated with CPI-7, sampled at days 0, 7, and 14."
    )

    if not _safe_exists(OVERVIEW_CSV) or not _safe_exists(BREAKDOWN_CSV):
        st.warning("Missing overview outputs. Run `make pipeline`.")
    else:
        summary = _read_csv(OVERVIEW_CSV).iloc[0]
        df_breakdown = _read_csv(BREAKDOWN_CSV)

        c1, c2, c3 = st.columns(3)
        c1.metric("Subjects", int(summary["total_subjects"]))
        c2.metric("Studies", int(summary["n_studies"]))
        c3.metric("Response rate", f"{summary['response_rate_pct']:.0f}%")

        st.divider()
        st.subheader("Headline finding")
        st.info("CD8 T cells at baseline are elevated in responders.")

        st.divider()
        st.subheader("Study breakdown")
        st.dataframe(
            df_breakdown.rename(columns={
                "condition": "Condition",
                "n_subjects": "Subjects",
                "n_responders": "Responders",
                "n_non_responders": "Non-responders",
            }),
            hide_index=True,
            width="stretch",
        )

        st.divider()
        with st.expander("How to use this dashboard"):
            st.markdown(
                """
- **Population Frequencies** — per-sample relative frequency (%) for each immune population.
- **Responder Analysis** — Mann–Whitney U + BH-corrected comparison of responders vs non-responders (CPI-7 / PBMC).
- **Baseline Cohort** — CPI-7 + PBMC subjects at day 0, with breakdowns by project, response, and sex.
- **Longitudinal Trajectories** — mean immune population percentages across days 0, 7, 14 by response group.
                """.strip()
            )

# ---------------------------------------------------------------------------
# Population Frequencies
# ---------------------------------------------------------------------------
with tab_frequencies:
    st.subheader("Frequency of each immune cell population in each sample")
    st.markdown(
        """
For each sample, total cell count is computed by summing counts across the five populations.
Relative frequency is reported as a percentage of total cells in that sample.
        """.strip()
    )

    if not _safe_exists(FREQ_CSV):
        st.warning("Missing `outputs/sample_population_frequencies.csv`. Run `python run_analysis.py`.")
    else:
        df_freq = _read_csv(FREQ_CSV)

        populations = sorted(df_freq["population"].unique().tolist())
        selected_pops = st.multiselect(
            "Population filter",
            options=populations,
            default=populations,
        )

        df_view = df_freq[df_freq["population"].isin(selected_pops)] if selected_pops else df_freq
        st.dataframe(df_view, hide_index=True, width="stretch")

# ---------------------------------------------------------------------------
# Responder Analysis
# ---------------------------------------------------------------------------
with tab_responder:
    st.subheader("Responder vs Non-responder Analysis (CPI-7 / PBMC)")

    st.markdown(
        """
Response status is determined after treatment based on clinical outcome criteria.
To explore whether immune composition patterns are associated with treatment response, we compare
relative frequencies of immune cell populations between eventual responders (yes) and non-responders (no)
across melanoma and NSCLC patients treated with **CPI-7**, using **PBMC** samples only.

**Statistical approach**
- Two-sided **Mann–Whitney U test** per population
- **Benjamini–Hochberg** correction across the five populations to control **false discovery rate (FDR)**
- Populations with **q < 0.05** are considered statistically significant
        """.strip()
    )

    time_label = st.selectbox(
        "Time from treatment start",
        options=["all", 0, 7, 14],
        index=0,
        key="responder_time",
    )
    label = str(time_label)

    df_resp = _load_responder_df(label)
    png_path = Path(str(RESPONDER_PNG_TMPL).format(label=label))

    if df_resp is None:
        st.warning(f"Missing responder results for time={label}. Run `python run_analysis.py`.")
    else:
        status, detail = _responder_headline(df_resp, label)
        st.info(status)
        if detail:
            st.caption(detail)

        if _safe_exists(png_path):
            st.image(str(png_path), use_column_width="always")
        else:
            st.warning(f"Missing plot for time={label}: {png_path.relative_to(BASE_DIR)}")

        with st.expander("Show statistical results table"):
            st.caption(
                "n_yes and n_no are the number of observations used in the test for each group "
                "at the selected timepoint. For pooled ('all') analysis, values represent per-subject "
                "means across available timepoints (each subject contributes one value per population)."
            )

            cols = [c for c in ["population", "n_yes", "n_no", "p_value", "q_value"] if c in df_resp.columns]
            show_df = df_resp[cols].copy()

            if "p_value" in show_df.columns:
                show_df["p_value"] = show_df["p_value"].map(lambda x: f"{float(x):.6f}")
            if "q_value" in show_df.columns:
                show_df["q_value"] = show_df["q_value"].map(lambda x: f"{float(x):.6f}")

            st.dataframe(show_df, hide_index=True, width="stretch")

# ---------------------------------------------------------------------------
# Baseline Cohort
# ---------------------------------------------------------------------------
with tab_baseline:
    st.subheader("Baseline Cohort (CPI-7 / PBMC / Day 0)")
    st.markdown(
        """
All subjects treated with **CPI-7** with **PBMC** samples at **baseline (time_from_treatment_start = 0)**,
across both melanoma and NSCLC conditions.
        """.strip()
    )

    if not _safe_exists(BASELINE_CSV):
        st.warning("Missing `outputs/baseline_cohort_summary.csv`. Run `python run_analysis.py`.")
    else:
        df_base = _read_csv(BASELINE_CSV)

        c1, c2 = st.columns(2)
        c1.metric("Baseline samples (PBMC)", len(df_base))
        c2.metric("Unique subjects", df_base["subject_id"].nunique())

        col_proj, col_resp, col_sex = st.columns(3)

        with col_proj:
            st.caption("Samples per project")
            proj_tbl = (
                df_base["project"]
                .fillna("unknown")
                .value_counts()
                .rename_axis("project")
                .reset_index(name="samples")
            )
            st.dataframe(proj_tbl, hide_index=True, width="stretch")

        with col_resp:
            st.caption("Subjects by response")
            resp_tbl = (
                df_base.drop_duplicates("subject_id")["response"]
                .fillna("unknown")
                .value_counts()
                .rename_axis("response")
                .reset_index(name="subjects")
            )
            st.dataframe(resp_tbl, hide_index=True, width="stretch")

        with col_sex:
            st.caption("Subjects by sex")
            sex_tbl = (
                df_base.drop_duplicates("subject_id")["sex"]
                .fillna("unknown")
                .value_counts()
                .rename_axis("sex")
                .reset_index(name="subjects")
            )
            st.dataframe(sex_tbl, hide_index=True, width="stretch")

        with st.expander("View full baseline sample list"):
            st.dataframe(df_base, hide_index=True, width="stretch")

# ---------------------------------------------------------------------------
# Longitudinal Trajectories
# ---------------------------------------------------------------------------
with tab_longitudinal:
    st.subheader("Longitudinal Immune Trajectories (CPI-7 / PBMC)")
    st.markdown(
        "Mean percentage of each immune population at days 0, 7, and 14, "
        "split by responder status across CPI-7-treated PBMC samples."
    )

    if not _safe_exists(LONGITUDINAL_CSV):
        st.warning("Missing `outputs/longitudinal_trajectories.csv`. Run `make pipeline`.")
    else:
        df_traj = _read_csv(LONGITUDINAL_CSV)

        condition_filter = st.selectbox(
            "Condition",
            options=["both", "melanoma", "nsclc"],
            index=0,
            key="longitudinal_condition",
        )

        if condition_filter == "both":
            df_plot = df_traj.groupby(
                ["timepoint", "response", "population"], as_index=False
            )["mean_percentage"].mean()
        else:
            df_plot = df_traj[df_traj["condition"] == condition_filter].copy()

        populations = sorted(df_plot["population"].unique())
        cols = st.columns(3)
        for i, pop in enumerate(populations):
            with cols[i % 3]:
                df_pop = df_plot[df_plot["population"] == pop]
                df_chart = (
                    df_pop.pivot(index="timepoint", columns="response", values="mean_percentage")
                    .sort_index()
                    .rename(columns={"yes": "Responders", "no": "Non-responders"})
                )
                df_chart.index.name = "Day"
                df_chart.columns.name = None
                st.caption(pop.replace("_", " ").title())
                st.line_chart(df_chart)