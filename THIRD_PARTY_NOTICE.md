# Third-party material included in this archive

## 1. Liu et al. (2025) supplementary figures

**File:** `core_analysis/evidence/Liu_2025_JCMM_70945_supplementary_figures.zip` (7.5 MB)

**Source:** Liu et al., *J Cell Mol Med*. 2025;29:e70945. DOI: 10.1111/jcmm.70945
Retrieved from Europe PMC supplementary files: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12611608/supplementaryFiles`

**Why it is here:** The supplementary figures S1–S4 were read directly to adjudicate the
TLN1/GYS1 claims (in particular the direction inconsistency between supplementary figure S3A and
the source paper's narrative of sustained up-regulation). This file is the provenance evidence
for those adjudication decisions.

**Action required before making this repository public:**

Confirm that redistribution of this file is permitted under the licence of the source article
(J Cell Mol Med; check whether the article is open access and under which licence). If
redistribution is not permitted — or if you prefer not to make the judgement call — remove the
file from the public repository and cite the source instead:

```bash
git rm --cached "core_analysis/evidence/Liu_2025_JCMM_70945_supplementary_figures.zip"
```

and replace it with a pointer file, for example
`core_analysis/evidence/Liu_2025_JCMM_70945_supplementary_figures.SOURCE.txt` containing:

```
Liu et al., J Cell Mol Med. 2025;29:e70945. DOI: 10.1111/jcmm.70945
Supplementary figures S1-S4 retrieved from Europe PMC:
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12611608/supplementaryFiles
Retrieved: <date>
Not redistributed here for licensing reasons.
```

The remainder of the archive does not depend on the binary file for reproducibility: all
adjudication decisions derived from it are already recorded in
`core_analysis/tables/source_adjudicated_subset_claim_disposition.csv` and
`core_analysis/evidence/external_evidence_decisions.csv`.

## 2. NCBI GEO metadata files

**Files:** `core_analysis/evidence/GSE*_metadata.soft`, `GSE54456_MAoverlappedsamples.txt.gz`

Public NCBI GEO records (GSE13355, GSE14905, GSE54456, GSE66511, GSE78097, GSE121212,
GSE201827, GSE295540). NCBI terms of use permit redistribution of retrieved records with
attribution; they are retained here as provenance evidence.

## 3. GSE295540 raw counts

**File:** `core_analysis/inputs/upload/GSE295540_raw_counts_All_samples.csv.gz` (1.0 MB)

Counts retrieved from the public GEO series, retained to make the RNA-seq boundary cohort
reproducible without a separate download step.
