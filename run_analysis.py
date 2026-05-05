import sqlite3
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "immune_trial.db"
OUTPUTS_DIR = BASE_DIR / "outputs"
RESPONDER_DIR = OUTPUTS_DIR / "responder_analysis"
POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")

# Allowed deviation from 100.0 when summing per-sample percentages (float rounding)
PCT_SUM_TOLERANCE = 0.01


def get_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def run_frequency_metrics(conn: sqlite3.Connection) -> None:
    """Generate outputs/sample_population_frequencies.csv with one row per (sample, population)."""
    query = """
        SELECT
            cc.sample_id AS sample,
            SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
            cc.population,
            cc.count
        FROM cell_counts cc
        ORDER BY cc.sample_id, cc.population
    """
    rows = conn.execute(query).fetchall()

    # --- Build output records and validate before writing ---
    records = []
    pct_sums: dict[str, float] = {}
    pop_counts: dict[str, int] = {}

    for row in rows:
        sample = row["sample"]
        total_count = row["total_count"]
        count = row["count"]

        if not total_count:
            raise RuntimeError(f"Zero total_count for sample {sample}")

        pct = count / total_count * 100.0
        records.append((sample, total_count, row["population"], count, round(pct, 6)))

        pct_sums[sample] = pct_sums.get(sample, 0.0) + pct
        pop_counts[sample] = pop_counts.get(sample, 0) + 1

    # Row count: must equal 5 * number_of_samples
    sample_count = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    expected_rows = sample_count * len(POPULATIONS)
    if len(records) != expected_rows:
        raise RuntimeError(
            f"Frequency metrics row count mismatch: got {len(records)}, expected {expected_rows} "
            f"(samples={sample_count}, populations={len(POPULATIONS)})"
        )

    # Each sample must have exactly 5 population rows
    bad_counts = {s: n for s, n in pop_counts.items() if n != len(POPULATIONS)}
    if bad_counts:
        raise RuntimeError(
            f"Samples with wrong population count (expected {len(POPULATIONS)}): {bad_counts}"
        )

    # Per-sample percentages must sum to ~100
    bad_sums = {s: v for s, v in pct_sums.items() if abs(v - 100.0) > PCT_SUM_TOLERANCE}
    if bad_sums:
        worst = max(bad_sums.items(), key=lambda x: abs(x[1] - 100.0))
        raise RuntimeError(
            f"Percentage sums deviate from 100 beyond tolerance={PCT_SUM_TOLERANCE} "
            f"in {len(bad_sums)} sample(s); worst: {worst[0]}={worst[1]:.6f}"
        )

    # --- All checks passed; write output ---
    out_path = OUTPUTS_DIR / "sample_population_frequencies.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "total_count", "population", "count", "percentage"])
        writer.writerows(records)

    print(f"[Frequency metrics] {out_path} written ({len(records)} rows, {sample_count} samples)")


def _bh_correction(p_values: list[float]) -> list[float]:
    """Benjamini–Hochberg FDR correction. Returns q-values in the same order as input."""
    n = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q_values = [0.0] * n
    running_min = 1.0
    for rev_rank, (orig_idx, p) in enumerate(reversed(indexed)):
        rank = n - rev_rank  # 1-based rank (descending order)
        q = p * n / rank
        running_min = min(running_min, q)
        q_values[orig_idx] = running_min
    return q_values


VALID_TIMEPOINTS = (0, 7, 14)

# ---------------------------------------------------------------------------
# Helpers for responder analysis
# ---------------------------------------------------------------------------

def _normalize_time_filter(time_filter: "int | str | None") -> "int | None":
    """Validate and coerce time_filter to an int (or None for pooled analysis).

    Accepted values:
        None, "all"          → returns None (no timepoint predicate)
        int or numeric str   → coerced to int, must be in VALID_TIMEPOINTS

    Raises ValueError for unrecognized or out-of-range values.
    """
    if time_filter is None or time_filter == "all":
        return None

    try:
        value = int(time_filter)
    except (TypeError, ValueError):
        raise ValueError(
            f"Invalid time_filter {repr(time_filter)}: expected None, 'all', "
            f"or an integer-compatible value."
        )

    if value not in VALID_TIMEPOINTS:
        raise ValueError(
            f"Invalid time_filter {repr(time_filter)}: {value} is not in "
            f"VALID_TIMEPOINTS {VALID_TIMEPOINTS}."
        )

    return value


def _fetch_cohort_rows(
    conn: sqlite3.Connection,
    time_filter: "int | str | None",
) -> list[sqlite3.Row]:
    """Return percentage rows for CPI-7+PBMC across all conditions, optionally by timepoint.

    Always includes subject_id so callers can aggregate to subject level if needed.
    time_filter is validated and normalized via _normalize_time_filter.
    """
    normalized = _normalize_time_filter(time_filter)
    base = """
        SELECT
            sub.subject_id,
            sub.response,
            cc.population,
            CAST(cc.count AS REAL) * 100.0 / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
        FROM samples sa
        JOIN subjects sub ON sa.subject_id = sub.subject_id
        JOIN cell_counts cc ON sa.sample_id = cc.sample_id
        WHERE sub.treatment = 'CPI-7'
          AND sa.sample_type = 'PBMC'
          AND sub.response IN ('yes', 'no')
    """
    order = " ORDER BY cc.population, sub.response"
    if normalized is None:
        return conn.execute(base + order).fetchall()
    return conn.execute(
        base + " AND sa.time_from_treatment_start = ?" + order,
        (normalized,),
    ).fetchall()


def _group_percentages(rows: list[sqlite3.Row]) -> dict[str, dict[str, list[float]]]:
    """Bucket sample-level percentage values by population → response group.

    Used for day-specific analyses where each row is one independent sample.
    """
    groups: dict[str, dict[str, list[float]]] = {
        pop: {"yes": [], "no": []} for pop in POPULATIONS
    }
    for row in rows:
        groups[row["population"]][row["response"]].append(row["percentage"])
    return groups


def _group_percentages_subject_mean(
    rows: list[sqlite3.Row],
) -> dict[str, dict[str, list[float]]]:
    """Aggregate to subject-level mean percentages by population → response group.

    Used for the pooled ("all") analysis to avoid treating repeated measures from
    the same subject as independent observations in the MWU test.

    Each subject contributes exactly one value per population (the mean across
    their available PBMC timepoints, typically 3: days 0, 7, 14).
    """
    from collections import defaultdict

    # (subject_id, population) → list of timepoint percentages
    subject_pcts: dict[tuple, list[float]] = defaultdict(list)
    # subject_id → response (validated for consistency)
    subject_response: dict[str, str] = {}

    for row in rows:
        sid = row["subject_id"]
        pop = row["population"]
        resp = row["response"]
        pct = row["percentage"]

        if sid in subject_response and subject_response[sid] != resp:
            raise RuntimeError(
                f"Subject '{sid}' appears with conflicting response values in cohort rows."
            )
        subject_response[sid] = resp
        subject_pcts[(sid, pop)].append(pct)

    # Validate: no subject contributes more than 3 timepoints per population
    over_limit = {
        k: v for k, v in subject_pcts.items() if len(v) > 3
    }
    if over_limit:
        raise RuntimeError(
            f"Subjects with >3 timepoints per population (unexpected): "
            f"{list(over_limit.keys())[:5]}"
        )

    groups: dict[str, dict[str, list[float]]] = {
        pop: {"yes": [], "no": []} for pop in POPULATIONS
    }
    for (sid, pop), pct_list in subject_pcts.items():
        resp = subject_response[sid]
        groups[pop][resp].append(sum(pct_list) / len(pct_list))

    return groups


def _run_mwu_with_bh(
    groups: dict[str, dict[str, list[float]]],
    label: str,
) -> list[dict]:
    """Mann–Whitney U + BH correction across all populations for one timepoint."""
    results = []
    for pop in POPULATIONS:
        yes_vals = groups[pop]["yes"]
        no_vals = groups[pop]["no"]
        n_yes, n_no = len(yes_vals), len(no_vals)
        if n_yes == 0 or n_no == 0:
            raise RuntimeError(
                f"[time={label}] Population '{pop}' has empty group: "
                f"n_yes={n_yes}, n_no={n_no}. Cannot run Mann–Whitney U test."
            )
        _, p_val = mannwhitneyu(yes_vals, no_vals, alternative="two-sided")
        results.append({"population": pop, "n_yes": n_yes, "n_no": n_no, "p_value": p_val})

    q_values = _bh_correction([r["p_value"] for r in results])
    for r, q in zip(results, q_values):
        r["q_value"] = q
    return results


def _write_differential_csv(results: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["population", "n_yes", "n_no", "p_value", "q_value"]
        )
        writer.writeheader()
        writer.writerows(results)


def _write_stratification_boxplot(
    groups: dict[str, dict[str, list[float]]],
    label: str,
    path: Path,
) -> None:
    time_str = "All Times from treatment start" if label == "all" else f"Time from treatment start = {label}"
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(4 * len(POPULATIONS), 5), sharey=False)
    for ax, pop in zip(axes, POPULATIONS):
        ax.boxplot(
            [groups[pop]["yes"], groups[pop]["no"]],
            labels=["Responder\n(yes)", "Non-responder\n(no)"],
            patch_artist=True,
            boxprops=dict(facecolor="#a8c8f0"),
            medianprops=dict(color="black", linewidth=2),
        )
        ax.set_title(pop.replace("_", " ").title(), fontsize=10)
        ax.set_ylabel("Percentage (%)")
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle(
        f"Responder vs Non-responder — CPI-7 / PBMC ({time_str})",
        fontsize=12,
        y=1.02,
    )
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main entry point for responder analysis (called once per timepoint)
# ---------------------------------------------------------------------------

def run_responder_analysis(
    conn: sqlite3.Connection,
    time_filter: "int | str | None" = None,
) -> None:
    """
    Responder vs non-responder differential analysis for CPI-7+PBMC across all conditions.

    time_filter:
        None / "all" — include all timepoints (pooled)
        0, 7, 14     — restrict to that samples.time_from_treatment_start value
    """
    label = "all" if (time_filter is None or time_filter == "all") else str(time_filter)

    rows = _fetch_cohort_rows(conn, time_filter)
    if not rows:
        raise RuntimeError(
            f"No rows for CPI-7+PBMC cohort (time={label!r}). "
            "Check data filters."
        )

    # Pooled analysis: average across timepoints per subject to preserve independence.
    # Day-specific analyses: each sample already represents one independent subject.
    if label == "all":
        groups = _group_percentages_subject_mean(rows)
    else:
        groups = _group_percentages(rows)
    results = _run_mwu_with_bh(groups, label)

    # Time-labelled outputs
    csv_path = RESPONDER_DIR / f"responder_differential_analysis_time_{label}.csv"
    _write_differential_csv(results, csv_path)
    print(f"[Responder analysis | time={label}] {csv_path} written")

    plot_path = RESPONDER_DIR / f"responder_stratification_boxplot_time_{label}.png"
    _write_stratification_boxplot(groups, label, plot_path)
    print(f"[Responder analysis | time={label}] {plot_path} saved")



def run_longitudinal_trajectories(conn: sqlite3.Connection) -> None:
    """Write outputs/longitudinal_trajectories.csv: mean % per population/timepoint/response/condition."""
    from collections import defaultdict

    rows = conn.execute("""
        SELECT
            sa.time_from_treatment_start AS timepoint,
            sub.response,
            sub.condition,
            cc.population,
            CAST(cc.count AS REAL) * 100.0 / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
        FROM samples sa
        JOIN subjects sub ON sa.subject_id = sub.subject_id
        JOIN cell_counts cc ON sa.sample_id = cc.sample_id
        WHERE sub.treatment = 'CPI-7'
          AND sa.sample_type = 'PBMC'
          AND sub.response IN ('yes', 'no')
        ORDER BY sa.time_from_treatment_start, sub.response, sub.condition, cc.population
    """).fetchall()

    if not rows:
        raise RuntimeError("No data found for longitudinal trajectories. Check CPI-7 / PBMC filters.")

    sums: dict[tuple, float] = defaultdict(float)
    counts: dict[tuple, int] = defaultdict(int)
    for row in rows:
        key = (row["timepoint"], row["response"], row["condition"], row["population"])
        sums[key] += row["percentage"]
        counts[key] += 1

    records = [
        {
            "timepoint": tp,
            "response": resp,
            "condition": cond,
            "population": pop,
            "mean_percentage": round(sums[(tp, resp, cond, pop)] / counts[(tp, resp, cond, pop)], 6),
        }
        for (tp, resp, cond, pop) in sorted(sums)
    ]

    out_path = OUTPUTS_DIR / "longitudinal_trajectories.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["timepoint", "response", "condition", "population", "mean_percentage"]
        )
        writer.writeheader()
        writer.writerows(records)

    print(f"[Longitudinal trajectories] {out_path} written ({len(records)} rows)")


def run_study_overview(conn: sqlite3.Connection) -> None:
    """Write outputs/study_overview.csv (summary) and outputs/study_breakdown.csv (per-condition)."""
    summary_row = conn.execute("""
        SELECT
            COUNT(*)                                                              AS total_subjects,
            COUNT(DISTINCT condition)                                             AS n_studies,
            SUM(CASE WHEN response = 'yes' THEN 1 ELSE 0 END)                    AS n_responders,
            ROUND(100.0 * SUM(CASE WHEN response = 'yes' THEN 1 ELSE 0 END) / COUNT(*), 1) AS response_rate_pct
        FROM subjects
        WHERE treatment = 'CPI-7'
    """).fetchone()

    summary_path = OUTPUTS_DIR / "study_overview.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["total_subjects", "n_studies", "n_responders", "response_rate_pct"])
        writer.writeheader()
        writer.writerow(dict(summary_row))
    print(f"[Study overview] {summary_path} written")

    breakdown_rows = conn.execute("""
        SELECT
            condition,
            COUNT(*)                                                 AS n_subjects,
            SUM(CASE WHEN response = 'yes' THEN 1 ELSE 0 END)       AS n_responders,
            SUM(CASE WHEN response = 'no'  THEN 1 ELSE 0 END)       AS n_non_responders
        FROM subjects
        WHERE treatment = 'CPI-7'
        GROUP BY condition
        ORDER BY condition
    """).fetchall()

    breakdown_path = OUTPUTS_DIR / "study_breakdown.csv"
    with breakdown_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["condition", "n_subjects", "n_responders", "n_non_responders"])
        writer.writeheader()
        writer.writerows([dict(r) for r in breakdown_rows])
    print(f"[Study overview] {breakdown_path} written")


def run_baseline_cohort_summary(conn: sqlite3.Connection) -> None:
    """Write outputs/baseline_cohort_summary.csv: all CPI-7 + PBMC subjects at day 0, both conditions."""
    rows = conn.execute("""
        SELECT
            sub.subject_id,
            sub.project,
            sub.condition,
            sub.sex,
            sub.response,
            sa.sample_id,
            sa.time_from_treatment_start,
            COALESCE(b_cell.count, 0)        AS b_cell_count,
            COALESCE(cd8_t_cell.count, 0)    AS cd8_t_cell_count,
            COALESCE(cd4_t_cell.count, 0)    AS cd4_t_cell_count,
            COALESCE(nk_cell.count, 0)       AS nk_cell_count,
            COALESCE(monocyte.count, 0)      AS monocyte_count
        FROM samples sa
        JOIN subjects sub ON sa.subject_id = sub.subject_id
        LEFT JOIN cell_counts b_cell     ON sa.sample_id = b_cell.sample_id     AND b_cell.population     = 'b_cell'
        LEFT JOIN cell_counts cd8_t_cell ON sa.sample_id = cd8_t_cell.sample_id AND cd8_t_cell.population = 'cd8_t_cell'
        LEFT JOIN cell_counts cd4_t_cell ON sa.sample_id = cd4_t_cell.sample_id AND cd4_t_cell.population = 'cd4_t_cell'
        LEFT JOIN cell_counts nk_cell    ON sa.sample_id = nk_cell.sample_id    AND nk_cell.population    = 'nk_cell'
        LEFT JOIN cell_counts monocyte   ON sa.sample_id = monocyte.sample_id   AND monocyte.population   = 'monocyte'
        WHERE sub.treatment = 'CPI-7'
          AND sa.sample_type = 'PBMC'
          AND sa.time_from_treatment_start = 0
        ORDER BY sub.condition, sub.subject_id
    """).fetchall()

    if not rows:
        raise RuntimeError("No baseline samples found for CPI-7 + PBMC at time=0.")

    out_path = OUTPUTS_DIR / "baseline_cohort_summary.csv"
    pop_cols = [f"{p}_count" for p in POPULATIONS]
    fieldnames = ["subject_id", "project", "condition", "sex", "response", "sample_id", "time_from_treatment_start"] + pop_cols
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fieldnames})

    by_condition: dict[str, int] = {}
    for row in rows:
        by_condition[row["condition"]] = by_condition.get(row["condition"], 0) + 1
    condition_summary = ", ".join(f"{c}: {n}" for c, n in sorted(by_condition.items()))
    print(f"[Baseline cohort summary] {out_path} written ({len(rows)} samples — {condition_summary})")


def main() -> None:
    OUTPUTS_DIR.mkdir(exist_ok=True)
    RESPONDER_DIR.mkdir(exist_ok=True)

    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found: {DB_PATH}. Run python load_data.py first."
        )

    with get_connection() as conn:
        run_frequency_metrics(conn)
        for time_filter in ("all", *VALID_TIMEPOINTS):
            run_responder_analysis(conn, time_filter)
        run_baseline_cohort_summary(conn)
        run_longitudinal_trajectories(conn)
        run_study_overview(conn)


if __name__ == "__main__":
    main()