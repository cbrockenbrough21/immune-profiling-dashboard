# Database Schema

`immune_trial.db` uses a three-table normalized schema stored in SQLite. The design separates subject-level metadata from sample-level metadata, and stores immune population measurements in long format to keep the schema stable as new cell types are added.

---

## Design Decisions

**Treatment is stored at the subject level.** It is a property of enrollment, not of any individual sample draw. Storing it on `subjects` avoids redundancy and makes cohort filtering by treatment straightforward without touching the samples table.

**Response is subject-level and nullable.** Clinical response is an outcome attributed to the patient as a whole, not to a specific sample or timepoint. It is stored on `subjects` and may be null for untreated subjects.

**Each subject belongs to one project.** Project membership is a fixed property of enrollment. The schema can be extended with a `projects` table and join table if subjects appear across multiple projects in future datasets.

**Cell counts are stored in long format.** Each (sample, population) pair is its own row rather than a separate column per population. Adding a new cell type requires no schema change, only new rows. It also makes population-level filtering and aggregation straightforward in SQL.

---

## subjects

One row per enrolled patient.

| Column | Type | Nullable | Description |
|--------|------|----------|-------------|
| subject_id | TEXT | NOT NULL (PK) | Unique patient identifier |
| project | TEXT | NOT NULL | Project identifier (e.g., PROJ-MEL-01) |
| condition | TEXT | NOT NULL | Clinical condition (e.g., melanoma, nsclc) |
| age | INTEGER | nullable | Age at enrollment |
| sex | TEXT | nullable | Biological sex (M/F) |
| treatment | TEXT | NOT NULL | Treatment assigned to the subject |
| response | TEXT | nullable | Clinical response outcome (yes/no); null for untreated subjects |

---

## samples

One row per biological sample. A subject with three timepoints has three rows here.

| Column | Type | Nullable | Description |
|--------|------|----------|-------------|
| sample_id | TEXT | NOT NULL (PK) | Unique sample identifier |
| subject_id | TEXT | NOT NULL (FK) | References subjects.subject_id |
| sample_type | TEXT | NOT NULL | Sample type (e.g., PBMC) |
| time_from_treatment_start | INTEGER | NOT NULL | Days from treatment start (0, 7, 14) |

`subject_id` is a foreign key into `subjects`. Cascading deletes are not enforced; the database is rebuilt from scratch on each pipeline run.

---

## cell_counts

One row per (sample, immune population) pair.

| Column | Type | Nullable | Description |
|--------|------|----------|-------------|
| sample_id | TEXT | NOT NULL (FK) | References samples.sample_id |
| population | TEXT | NOT NULL | Immune cell type (b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte) |
| count | INTEGER | NOT NULL | Raw cell count |

Primary key: `(sample_id, population)` — enforces one count per population per sample and makes point lookups efficient.

Relative frequencies are computed at query time rather than stored. This keeps the source of truth in raw counts and avoids stale derived values if counts change.

---

## Indexes

| Index | Table | Columns | Purpose |
|-------|-------|---------|---------|
| idx_subjects_project | subjects | project | Cohort filtering by project |
| idx_subjects_treatment | subjects | treatment | Cohort filtering by treatment |
| idx_samples_subject | samples | subject_id | Join from subjects to samples |
| idx_samples_type_time | samples | sample_type, time_from_treatment_start | Filtering to PBMC samples at specific timepoints |
| idx_cell_counts_population | cell_counts | population | Filtering to specific immune populations |

Indexes cover the primary join and filter paths used in the analytical queries. The `(sample_type, time_from_treatment_start)` composite index is particularly useful for the baseline and timepoint-stratified cohort queries that always apply both predicates together.
