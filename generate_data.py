"""
Generates synthetic PBMC immune profiling data for two clinical trial cohorts.

One embedded signal: CD8 T cells are meaningfully elevated in responders at
baseline (day 0). All other populations and all non-baseline timepoints are noise.

Output: data/cell-count.csv
"""
import csv
import random
from pathlib import Path

SEED = 42
OUTPUT_PATH = Path("data/cell-count.csv")

N_SUBJECTS = 80
RESPONDER_RATE = 0.40
TIMEPOINTS = [0, 7, 14]
POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")
TREATMENT = "CPI-7"
SAMPLE_TYPE = "PBMC"

STUDIES = [
    {"project": "PROJ-MEL-01", "condition": "melanoma", "prefix": "MEL", "n": 40},
    {"project": "PROJ-NSC-01", "condition": "nsclc",    "prefix": "NSC", "n": 40},
]

# Noise distributions (mean, std) — identical for responders and non-responders
NOISE = {
    "b_cell":     (8_000,  2_000),
    "cd4_t_cell": (25_000, 5_000),
    "nk_cell":    (8_000,  2_000),
    "monocyte":   (15_000, 3_000),
}

# CD8 signal at day 0 only; all other timepoints use the noise distribution
CD8_RESPONDER_BASELINE    = (18_000, 3_500)
CD8_NONRESPONDER_BASELINE = (12_000, 3_500)
CD8_NOISE                 = (15_000, 5_000)


def pos_int(mean: float, std: float) -> int:
    return max(1, round(random.gauss(mean, std)))


def subject_rows(subject_id: str, project: str, condition: str, is_responder: bool) -> list[dict]:
    age = random.randint(45, 75)
    sex = random.choice(["M", "F"])
    response = "yes" if is_responder else "no"

    rows = []
    for day in TIMEPOINTS:
        counts = {}
        for pop in POPULATIONS:
            if pop == "cd8_t_cell":
                if day == 0:
                    params = CD8_RESPONDER_BASELINE if is_responder else CD8_NONRESPONDER_BASELINE
                else:
                    params = CD8_NOISE
            else:
                params = NOISE[pop]
            counts[pop] = pos_int(*params)

        rows.append({
            "subject":                   subject_id,
            "sample":                    f"{subject_id}-D{day:02d}",
            "project":                   project,
            "condition":                 condition,
            "treatment":                 TREATMENT,
            "sample_type":               SAMPLE_TYPE,
            "time_from_treatment_start": day,
            "response":                  response,
            "age":                       age,
            "sex":                       sex,
            **counts,
        })
    return rows


def main() -> None:
    random.seed(SEED)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    n_responders_per_study = round(N_SUBJECTS * RESPONDER_RATE) // len(STUDIES)  # 16

    all_rows: list[dict] = []
    for study in STUDIES:
        n = study["n"]
        flags = [True] * n_responders_per_study + [False] * (n - n_responders_per_study)
        random.shuffle(flags)

        for i, is_responder in enumerate(flags):
            subject_id = f"{study['prefix']}-{i + 1:03d}"
            all_rows.extend(
                subject_rows(subject_id, study["project"], study["condition"], is_responder)
            )

    fieldnames = [
        "subject", "sample", "project", "condition", "treatment",
        "sample_type", "time_from_treatment_start", "response",
        "age", "sex",
        *POPULATIONS,
    ]

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    n_responders = sum(1 for r in all_rows if r["response"] == "yes" and r["time_from_treatment_start"] == 0)
    n_nonresponders = sum(1 for r in all_rows if r["response"] == "no" and r["time_from_treatment_start"] == 0)
    print(f"Generated {len(all_rows)} rows ({N_SUBJECTS} subjects × {len(TIMEPOINTS)} timepoints)")
    print(f"  Responders: {n_responders}  |  Non-responders: {n_nonresponders}")
    print(f"  Output: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
