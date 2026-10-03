# Psoriasis transcriptomic reproducibility archive v1.0

This archive contains the processed inputs, evidence files, analysis code, numerical outputs, sensitivity analyses, paired audit-impact analyses, and publication figures supporting the study:

**Dataset provenance, sample overlap and direction preservation in cross-cohort reproducibility of transcriptomic biomarker claims: an empirical methodological study in psoriasis**

## What this archive represents

The final source-adjudicated development analysis contains **33 paper-gene claims representing 32 genes and 176 eligible cohort contributions**. The archived inventory contains **39 claims representing 36 genes**. A separate RNA-seq cohort (GSE295540; 7 psoriasis and 6 healthy samples) is retained as a small external boundary case.

The current analysis has four main audit components:

1. **Source-use/provenance audit**: source-used cohorts are excluded on a paper-gene basis.
2. **Cross-accession sample-overlap correction**: the 42 samples shared between GSE13355 and GSE54456 are removed from every eligible GSE54456 contribution in the final synthesis, leaving 72 psoriasis and 60 healthy samples.
3. **Direction freezing**: published expression direction is fixed before benchmark interpretation; post hoc `max(AUC, 1-AUC)` reversal is not used in the principal analysis.
4. **Direction-preserving synthesis and uncertainty analysis**: claim-specific random-effects synthesis, modified Hartung-Knapp intervals, prediction intervals, leave-one-cohort-out checks, platform sensitivity, and related analyses.

A later paired audit-impact analysis quantifies what each of these methodological choices changed. It is explicitly post hoc and descriptive.

## Directory structure

- `core_analysis/code/` — executable Python analysis and verification scripts.
- `core_analysis/inputs/` — processed development matrices, frozen analysis inputs, overlap mappings/exclusion sets, and GSE295540 counts/phenotype files.
- `core_analysis/evidence/` — GEO metadata, official overlap evidence, source-paper supplementary evidence, preprocessing checks, and source-decision records.
- `core_analysis/results/` — numerical results from the verified core analysis release.
- `core_analysis/tables/` — detailed scientific analysis tables and adjudicated cohort/claim outputs.
- `core_analysis/figures/` — publication figures in PDF/PNG/SVG, including `FigureS1_Adjudicated_Forest.*`.
- `core_analysis/audit_impact/` — paired provenance, overlap, direction-freezing, and inference-layer outputs added on 2026-09-17.
- `core_analysis/legacy_patches/` — retained audit trail for earlier numerical-function repairs; not the primary execution path.
- `historical_notes/` — historical manuscript snapshots, manuscript-generation helpers, and older explanatory notes retained for provenance only. These are **not** the authoritative current journal manuscript.

## Primary result files

For the current manuscript, the authoritative development population is the **source-adjudicated 33-claim set**. In particular:

- `core_analysis/tables/source_adjudicated_subset_claim_disposition.csv`
- `core_analysis/tables/source_adjudicated_subset_cohort_inputs.csv`
- `core_analysis/tables/source_adjudicated_subset_development_meta.csv`
- `core_analysis/tables/source_adjudicated_subset_postholdout_meta.csv`

The files prefixed `confirmed_use_corrections_` and the original 39-claim grid are retained for audit/history and should not be substituted for the current principal development population.

Paired audit-impact outputs are in:

- `core_analysis/audit_impact/source_use_paired_33_claims.csv`
- `core_analysis/audit_impact/provenance_not_sufficiently_evaluable.csv`
- `core_analysis/audit_impact/overlap_correction_paired_30_claims.csv`
- `core_analysis/audit_impact/direction_freezing_counterfactual_176_contributions.csv`
- `core_analysis/audit_impact/inference_layer_paired_33_claims.csv`
- `core_analysis/audit_impact/audit_impact_summary.json`
- `core_analysis/audit_impact/main_table4_audit_impact.csv`

## Reproducing the analysis

The verified environment used Python 3.12 with dependencies listed in `core_analysis/requirements.txt`.

From the archive root, either run the two analysis layers separately:

```bash
cd core_analysis
python -m venv .venv
# activate the environment using the command appropriate for your operating system
pip install -r requirements.txt
python code/run_pipeline.py
python code/audit_impact_analysis.py
```

or, after installing the same requirements, run:

```bash
python run_all.py
```

`run_pipeline.py` verifies frozen input hashes before running the core computation. `audit_impact_analysis.py` regenerates the paired methodological comparisons in `core_analysis/audit_impact/`.

### Important note on historical manuscript helper files

The verified 2026-09-08 computational release originally included manuscript-generation helper code and manuscript-facing CSV tables. The journal manuscript was subsequently revised for BMC Medical Research Methodology, including updated overlap wording and the new audit-impact analysis. To avoid presenting outdated journal-facing text as current, those manuscript snapshots/helpers are intentionally omitted from this public reproducibility archive. Historical explanatory notes and the earlier release manifest are retained in `historical_notes/` for provenance.

## Key current quantitative checks

- Source-use adjudication removed **17 source-used cohort tests**.
- Five PANoptosis claims were left with only one eligible development cohort and are classified as **not sufficiently evaluable**, not as replication failures.
- All final eligible GSE54456 contributions use the overlap-corrected sample set: **72 psoriasis + 60 healthy** after removal of **42 shared samples**.
- The paired overlap analysis affects 30 claims but changes no normal-CI, mHK-CI, normal-PI, or t-based-PI classification.
- Five of 176 final cohort contributions have directional AUC below 0.5; all are TLN1. The `max(AUC,1-AUC)` analysis is included only as an illustrative counterfactual.
- Final development inference: conventional normal CI supports 32 positive claims and one negative claim; mHK supports 30 positive claims with CD28, ITGAL, and TLN1 crossing zero.

## Data provenance

Original expression datasets are publicly available from NCBI GEO under:

GSE13355, GSE14905, GSE54456, GSE66511, GSE78097, GSE121212, GSE201827, and GSE295540.

This archive contains processed analysis matrices and metadata/evidence needed for the reported claim-level analyses. It does not repackage the full original GEO repositories or claim to reprocess all raw CEL/FASTQ data.

## Integrity

`SHA256SUMS.txt` at the archive root lists SHA-256 hashes for every file included in this release (excluding the checksum file itself).

## Versioning

Archive version: **1.0**  
Prepared: **2026-09-17**

If analysis files are changed after public release, create a new repository version rather than silently replacing this archive.

## License

No new license is assigned by this archive-generation step. The authors should select an appropriate license when publishing the record in Zenodo and ensure that any redistributed third-party material remains compatible with its original license and attribution requirements.
