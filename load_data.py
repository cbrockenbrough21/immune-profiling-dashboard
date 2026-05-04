import csv
import sqlite3
from pathlib import Path
from typing import Optional


DB_FILENAME = "immune_trial.db"
CSV_CANDIDATES = ("cell-count.csv", "data/cell-count.csv")

# Canonical populations for this dataset (allowlist)
POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")

# Subject-level fields that must be consistent across all rows for the same subject
SUBJECT_FIELDS = ("project", "condition", "treatment", "age", "sex")


def normalize_optional_str(value: Optional[str]) -> Optional[str]:
    """
    Normalize CSV fields where missing values may appear as empty strings or whitespace.
    Returns None for empty/whitespace values.
    """
    if value is None:
        return None
    value = value.strip()
    return value or None


def parse_optional_int(value: Optional[str]) -> Optional[int]:
    """
    Parse an optional integer field from CSV. Returns None if missing/blank.
    """
    v = normalize_optional_str(value)
    if v is None:
        return None
    return int(v)


def resolve_csv_path(repo_root: Path) -> Path:
    for candidate in CSV_CANDIDATES:
        candidate_path = repo_root / candidate
        if candidate_path.exists():
            return candidate_path
    raise FileNotFoundError(
        f"Could not find input CSV. Checked: {', '.join(CSV_CANDIDATES)}"
    )


def rebuild_database(db_path: Path) -> sqlite3.Connection:
    if db_path.exists():
        db_path.unlink()
    connection = sqlite3.connect(db_path)
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def create_schema(connection: sqlite3.Connection) -> None:
    # NOTE: response is nullable (TEXT) to support untreated subjects (e.g., treatment='none')
    connection.executescript(
        """
        CREATE TABLE subjects (
            subject_id TEXT PRIMARY KEY NOT NULL,
            project TEXT NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER,
            sex TEXT,
            treatment TEXT NOT NULL,
            response TEXT
        );

        CREATE TABLE samples (
            sample_id TEXT PRIMARY KEY NOT NULL,
            subject_id TEXT NOT NULL,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL,
            FOREIGN KEY(subject_id) REFERENCES subjects(subject_id)
        );

        CREATE TABLE cell_counts (
            sample_id TEXT NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            PRIMARY KEY (sample_id, population),
            FOREIGN KEY(sample_id) REFERENCES samples(sample_id)
        );

        CREATE INDEX idx_subjects_project ON subjects(project);
        CREATE INDEX idx_subjects_treatment ON subjects(treatment);
        CREATE INDEX idx_samples_subject ON samples(subject_id);
        CREATE INDEX idx_samples_type_time ON samples(sample_type, time_from_treatment_start);
        CREATE INDEX idx_cell_counts_population ON cell_counts(population);
        """
    )


def ingest_csv(connection: sqlite3.Connection, csv_path: Path) -> tuple[int, int, int]:
    """
    Loads data from CSV into subjects, samples, and cell_counts.

    Response consistency rule:
    - Normalize blanks to NULL
    - Only enforce consistency across NON-NULL responses per subject
    - If a subject has multiple distinct non-null responses -> error
    """
    # Track a single canonical non-null response per subject (if present)
    subject_to_response: dict[str, Optional[str]] = {}

    # Track subject-level metadata for consistency validation
    subject_metadata: dict[str, dict[str, str]] = {}

    # Track sample_ids to detect duplicates in the CSV
    seen_sample_ids: set[str] = set()

    # subjects table rows (dedup by subject_id)
    subjects: dict[str, tuple[str, str, str, Optional[int], Optional[str], str, Optional[str]]] = {}

    # samples table rows (one per CSV row)
    samples: list[tuple[str, str, str, int]] = []

    # cell_counts table rows (one per population per sample)
    cell_counts: list[tuple[str, str, int]] = []

    with csv_path.open("r", newline="", encoding="utf-8-sig") as csv_file:
        reader = csv.DictReader(csv_file)

        # Minimal header validation for expected columns
        required_cols = {"subject", "sample", "project", "condition", "treatment", "sample_type", "time_from_treatment_start"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV is missing required columns: {sorted(missing)}")

        for pop in POPULATIONS:
            if reader.fieldnames is None or pop not in reader.fieldnames:
                raise ValueError(f"CSV is missing expected population column: {pop}")

        for row in reader:
            subject_id = row["subject"].strip()
            sample_id = row["sample"].strip()

            # Validate no duplicate sample_ids in the CSV
            if sample_id in seen_sample_ids:
                raise ValueError(f"Duplicate sample_id in CSV: {sample_id}")
            seen_sample_ids.add(sample_id)

            # Normalize optional fields
            response = normalize_optional_str(row.get("response"))
            sex = normalize_optional_str(row.get("sex"))
            age = parse_optional_int(row.get("age"))

            # Only compare non-null responses — nulls are allowed (e.g. untreated subjects)
            existing_response = subject_to_response.get(subject_id)

            if existing_response is None:
                # We haven't recorded a non-null response yet (or haven't seen this subject)
                if response is not None:
                    subject_to_response[subject_id] = response
                else:
                    # Keep as None for now
                    subject_to_response.setdefault(subject_id, None)
            else:
                # We already have a recorded value (could be None or a non-null response)
                if response is not None:
                    if existing_response is None:
                        subject_to_response[subject_id] = response
                    elif existing_response != response:
                        raise ValueError(
                            "Inconsistent response for subject "
                            f"{subject_id}: {existing_response} vs {response}"
                        )

            # Validate subject-level metadata consistency across rows
            current_meta = {
                "project": row["project"].strip(),
                "condition": row["condition"].strip(),
                "treatment": row["treatment"].strip(),
                "age": row.get("age", "").strip(),
                "sex": row.get("sex", "").strip(),
            }
            if subject_id not in subject_metadata:
                subject_metadata[subject_id] = current_meta
            else:
                stored = subject_metadata[subject_id]
                for field in SUBJECT_FIELDS:
                    if stored[field] != current_meta[field]:
                        raise ValueError(
                            f"Inconsistent {field} for subject {subject_id}: "
                            f"{stored[field]!r} vs {current_meta[field]!r}"
                        )

            # Insert subject once
            if subject_id not in subjects:
                subjects[subject_id] = (
                    subject_id,
                    row["project"].strip(),
                    row["condition"].strip(),
                    age,
                    sex,
                    row["treatment"].strip(),
                    subject_to_response.get(subject_id),
                )
            else:
                # If we later learned a non-null response, update the stored subject row
                # (ensures subjects.response reflects the canonical non-null response if present)
                current = subjects[subject_id]
                canonical_response = subject_to_response.get(subject_id)
                if current[-1] != canonical_response:
                    subjects[subject_id] = (*current[:-1], canonical_response)

            # Sample row
            samples.append(
                (
                    sample_id,
                    subject_id,
                    row["sample_type"].strip(),
                    int(row["time_from_treatment_start"]),
                )
            )

            # Cell counts
            for population in POPULATIONS:
                cell_counts.append((sample_id, population, int(row[population])))

    connection.executemany(
        """
        INSERT INTO subjects (
            subject_id, project, condition, age, sex, treatment, response
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        list(subjects.values()),
    )

    connection.executemany(
        """
        INSERT INTO samples (
            sample_id, subject_id, sample_type, time_from_treatment_start
        ) VALUES (?, ?, ?, ?)
        """,
        samples,
    )

    connection.executemany(
        """
        INSERT INTO cell_counts (sample_id, population, count)
        VALUES (?, ?, ?)
        """,
        cell_counts,
    )

    return len(subjects), len(samples), len(cell_counts)


def fetch_table_count(connection: sqlite3.Connection, table_name: str) -> int:
    return connection.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]


def validate_integrity(connection: sqlite3.Connection) -> None:
    sample_rows = fetch_table_count(connection, "samples")
    cell_count_rows = fetch_table_count(connection, "cell_counts")
    expected_cell_count_rows = sample_rows * len(POPULATIONS)

    if cell_count_rows != expected_cell_count_rows:
        raise ValueError(
            "Integrity check failed: total cell_counts rows must equal "
            f"samples * populations ({sample_rows} * {len(POPULATIONS)} = "
            f"{expected_cell_count_rows}), got {cell_count_rows}"
        )

    invalid_sample = connection.execute(
        """
        SELECT sample_id, COUNT(*) AS population_count
        FROM cell_counts
        GROUP BY sample_id
        HAVING population_count != ?
        LIMIT 1
        """,
        (len(POPULATIONS),),
    ).fetchone()
    if invalid_sample is not None:
        sample_id, population_count = invalid_sample
        raise ValueError(
            "Integrity check failed: each sample must have exactly "
            f"{len(POPULATIONS)} populations; sample {sample_id} has "
            f"{population_count}"
        )

    placeholders = ", ".join("?" for _ in POPULATIONS)
    invalid_population = connection.execute(
        f"""
        SELECT population
        FROM cell_counts
        WHERE population NOT IN ({placeholders})
        LIMIT 1
        """,
        POPULATIONS,
    ).fetchone()
    if invalid_population is not None:
        raise ValueError(f"Unexpected population value: {invalid_population[0]}")


def run() -> None:
    repo_root = Path(__file__).resolve().parent
    db_path = repo_root / DB_FILENAME
    csv_path = resolve_csv_path(repo_root)

    connection = rebuild_database(db_path)
    try:
        with connection:
            create_schema(connection)
            subject_count, sample_count, cell_count = ingest_csv(connection, csv_path)
            validate_integrity(connection)
    finally:
        connection.close()

    print("load_data.py completed successfully")
    print(f"Database: {db_path}")
    print(f"Input CSV: {csv_path}")
    print(f"subjects inserted: {subject_count}")
    print(f"samples inserted: {sample_count}")
    print(f"cell_counts inserted: {cell_count}")
    print("Integrity checks passed")


if __name__ == "__main__":
    run()