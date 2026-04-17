# Repeat Expansion Analysis Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Run with Singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![SLURM](https://img.shields.io/badge/executor-SLURM-blue.svg)](https://slurm.schedmd.com/)


A production-grade, HPC-scale Nextflow DSL2 pipeline for detecting, annotating and reporting short tandem repeat (STR) expansions in whole-exome sequencing (WES) cohorts.

Developed for the **RCSI epilepsy cohort** at RCSI designed to be reusable for any WES or WGS cohort with PED-defined family structure.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Files](#input-files)
- [Parameters](#parameters)
- [Output Structure](#output-structure)
- [Interpreting the Excel Report](#interpreting-the-excel-report)
- [REViewer : On-Demand Visualisation](#reviewer--on-demand-visualisation)
- [Resource Configuration](#resource-configuration)
- [Limitations](#limitations)
- [Citation](#citation)

---

## Overview

Standard SNV/indel variant calling pipelines do not detect short tandem repeat expansions. This pipeline fills that gap by:

1. **Extracting** reads overlapping 66 STR loci from full WES BAM files (hybrid seeking mode)
2. **Genotyping** repeat expansions using [ExpansionHunter v5](https://github.com/Illumina/ExpansionHunter)
3. **Annotating** per-sample VCFs with disease metadata via the [STRipy](https://stripy.org) API
4. **Scoring** each allele with population Z-scores, outlier flags, and allele direction
5. **Assessing de novo** status for all trio and duo families
6. **Reporting** one colour-coded Excel workbook per family  ready for clinical review

> **Important distinction:** The `Allele_DeNovo` flag indicates a repeat was not observed in either parent (inheritance finding). Whether a repeat is disease-causing is determined by the `PathogenicCutoff` column  a de novo repeat below the pathogenic threshold is not automatically pathogenic and must be interpreted in full clinical context.

---

## Pipeline Architecture

```
WES BAM files 
        │
        ▼
SAMTOOLS_EXTRACT_REGIONS     ← samtools view -L repeat_regions.bed
        │  mini-BAM 
        ▼
EXPANSIONHUNTER              ← v5.0.0, seeking mode, 1 CPU per sample
        │  per-sample: .vcf  .json  _realigned.sorted.bam
        ▼
SAMTOOLS_INDEX_REALIGNED     ← sort + index realigned BAM for REViewer
        │
        ▼
STRIPY_FETCH_LOCUS_REF       ← bootstrap once from /locus API → JSON cache
        │
        ▼
STRIPY_ANNOTATE_VCF          ← POST to /annotateVCF API per sample
        │  {sample}.Reannotated.vcf
        ▼
STRIPY_BUILD_SAMPLE_REPORT   ← /compare API + HPO lookup + de novo logic
        │  {sample}.stripy_report.tsv  (48 columns)
        ▼
STRIPY_AGGREGATE_COHORT_REPORT  ← one Excel per family
           {FamilyID}_{ProbandID}.xlsx
           Sheets: Summary | Proband | Father | Mother | Legend
```

### Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| **Hybrid seeking mode** | Region extraction → mini-BAM → EH seeking. Same speed as streaming but produces realigned BAMs for REViewer visualisation |
| **Static vs dynamic API** | STRipy `/locus` data fetched once (bootstrap). `/compare` called per sample only  eliminates ~16,000 redundant API calls across the cohort |
| **Per-allele columns** | All Z-score, outlier, direction and de novo fields are explicit A1/A2 columns  enables direct Excel filtering without text parsing |
| **OMIM as HPO bridge** | STRipy disease IDs (`HMNR7` etc.) are internal codes. OMIM IDs from the `/locus` reference bridge to `genes_to_phenotype.txt` for HPO lookup |
| **Separate containers** | `expansionhunter.sif` does not include samtools  sort/index runs in `samtools_1.23.sif` |

---

## Requirements

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| [Nextflow](https://nextflow.io) | ≥ 25.04.0 | Pipeline orchestration |
| [Singularity](https://sylabs.io) | ≥ 3.8 | Container runtime |
| [Java](https://openjdk.org) | 17 | Nextflow runtime |
| SLURM | any | HPC job scheduler |

### Singularity Images

Place all images in a shared directory (default: `/home/data/human_genetics/ONE4ALL/Docker_Images/`):

```
samtools_1.23.sif
expansionhunter.sif          # ExpansionHunter v5.0.0
python_3.14.2.sif            # with openpyxl installed
```

> **Note:** The Python container must have `openpyxl` installed. If your base Python image lacks it, build a derived image:
> ```bash
> cat > python_openpyxl.def << 'EOF'
> Bootstrap: localimage
> From: /path/to/python_3.14.2.sif
> %post
>     pip install openpyxl --break-system-packages
> EOF
> singularity build python_openpyxl.sif python_openpyxl.def
> ```

### Reference Files

| File | Notes |
|------|-------|
| `hg38.fasta` + `.fai` | Reference genome |
| `variant_catalog_v2_extended_offtarget_added_chr.json` | ExpansionHunter variant catalog (66 loci) |
| `repeat_regions.bed` | Padded BED of repeat loci  generate once (see below) |
| `hpo_genes_to_phenotype.txt` | [HPO](https://hpo.jax.org/data/annotations) gene-to-phenotype annotations |
| `stripy_locus_reference.json` | STRipy /locus cache  auto-generated on first pipeline run |

---

## Installation

```bash
# Clone the repository
git clone https://github.com/your-org/lighthouse-str-pipeline.git
cd lighthouse-str-pipeline

# Generate the repeat regions BED file (run once)
python3 scripts/generate_regions_bed.py \
    --catalog  /path/to/variant_catalog_v2_extended_offtarget_added_chr.json \
    --output   /path/to/reference/repeat_regions.bed \
    --padding  1000

```

The `stripy_locus_reference.json` bootstrap cache will be generated automatically on the first pipeline run. It fetches metadata for all 66 loci from the STRipy `/locus` API and saves permanently to `params.stripy_locus_ref`. Subsequent runs use the cached file  delete it to force a refresh.

---

## Quick Start

```bash
# Edit the SLURM launcher with your paths
vim run_expansionhunter.slurm

# Submit to SLURM
sbatch run_expansionhunter.slurm
```

Or run directly with Nextflow:

```bash
nextflow run main.nf \
    -profile slurm \
    -resume \
    --ped             /path/to/cohort.ped \
    --bam_dir         /path/to/bam/markdup \
    --fasta           /path/to/hg38.fasta \
    --variant_catalog /path/to/variant_catalog_v2.json \
    --regions_bed     /path/to/repeat_regions.bed \
    --outdir          /path/to/results \
    --stripy_locus_ref /path/to/stripy_locus_reference.json \
    --hpo_file        /path/to/hpo_genes_to_phenotype.txt \
 ```

Resume after failure  Nextflow caches all completed tasks:

```bash
nextflow run main.nf -profile slurm -resume [all other params]
```

---

## Input Files

### PED File

Standard 6-column PED format:

```
FamilyID  SampleID       FatherID       MotherID       Sex  Phenotype
LH0001    LH0001A_S4     LH0001B_S5     LH0001C_S2     1    2
LH0001    LH0001B_S5     0              0              1    0
LH0001    LH0001C_S2     0              0              2    0
```

- `Sex`: 1 = male, 2 = female
- `Phenotype`: 2 = affected (proband), 0/1 = unaffected
- Proband identified by `phenotype = 2`. Fallback: suffix-based detection (A/P = proband, B/D/F = father, C/M = mother)

### BAM Files

- Coordinate-sorted, duplicate-marked BAM files with index (`.bai`)
- Default suffix: `.md.bam` (override with `--bam_suffix`)
- Files named `{SampleID}{bam_suffix}`  must match PED column 2 exactly

---

## Parameters

### Core

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ped` | `null` | PED file (required) |
| `--bam_dir` | `null` | Directory containing BAM files (required) |
| `--fasta` | `null` | Reference FASTA (required) |
| `--variant_catalog` | `null` | EH variant catalog JSON (required) |
| `--regions_bed` | `null` | Padded repeat loci BED (required for seeking mode) |
| `--outdir` | `./results` | Output directory |
| `--bam_suffix` | `.md.bam` | BAM filename suffix |
| `--analysis_mode` | `seeking` | `seeking` (recommended) or `streaming` |

### STRipy Post-processing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--stripy_locus_ref` | *(path)* | STRipy locus reference JSON  auto-generated on first run |
| `--hpo_file` | *(path)* | HPO genes_to_phenotype.txt |
| `--cohort_id` | `LH_WES` | Cohort identifier used in output filenames |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `30` | Hard CPU cap across all jobs |
| `--max_memory` | `128.GB` | Hard memory cap |
| `--slurm_queue_size` | `25` | Max concurrent SLURM jobs |
| `--slurm_partition` | `ProdQ` | SLURM partition |

Override any parameter at the command line with `--parameter value`.

---

## Output Structure

```
results/
├── expansionhunter/
│   └── {SampleID}/
│       ├── {SampleID}.vcf                       # EH genotypes
│       ├── {SampleID}.json                      # EH full output
│       ├── {SampleID}_realigned.bam             # Unsorted realigned BAM
│       ├── {SampleID}_realigned.sorted.bam      # Sorted  use for REViewer
│       ├── {SampleID}_realigned.sorted.bam.bai  # BAM index
│       └── versions.yml
│
├── stripy/
│   ├── annotated_vcf/
│   │   └── {SampleID}.Reannotated.vcf           # STRipy-annotated VCF
│   └── reports/
│       └── {SampleID}.stripy_report.tsv         # 48-column per-sample TSV
│
├── cohort/
│   └── {FamilyID}_{ProbandID}.xlsx              # Per-family Excel workbook
│
└── pipeline_info/
    ├── completed_samples.tsv                    # Run completion manifest
    ├── timeline_*.html                          # Nextflow timeline
    ├── report_*.html                            # Nextflow execution report
    ├── trace_*.txt                              # Per-task trace
    └── dag_*.svg                                # Pipeline DAG
```

### 48-Column TSV Schema

Each `.stripy_report.tsv` contains one row per disease association per locus:

**Section A  Locus identity & clinical context (30 columns)**
`SampleID · Chr · Start · End · LocationCoordinates · Ref · Alt · GT · LocusCov · RepeatLength · FilterStatus · CallStatus · CallQuality · Locus · LocationRegion · RepeatType · Motif · Normal Range (Min/Max) · IntermediateRange (Min/Max) · PathogenicCutoff · Disease Name · Disease OMIM · Inheritance · Onset · Stripy_Gene · HPO_Gene · HPO Phenotype · HPO_Filtered`

**Section B  Per-allele population statistics (12 columns)**
`Allele1_Repeats · Allele2_Repeats · Literature:DiseaseName · Literature:Inheritance · Allele1_Range · Allele2_Range · Allele1_Zscore · Allele2_Zscore · Allele1_Outlier · Allele2_Outlier · Allele1_Direction · Allele2_Direction`

**Section C  Per-allele de novo assessment (6 columns)**
`Allele1_DeNovo · Allele1_DeNovo_Direction · Allele2_DeNovo · Allele2_DeNovo_Direction · DeNovo_Status · Caller`

---

## Interpreting the Excel Report

Each family receives one `.xlsx` workbook with five sheets:

| Sheet | Contents | Start here? |
|-------|----------|-------------|
| **Summary** | Proband rows flagged as outliers or de novo candidates only | ★ Yes |
| **Proband** | All 66 loci for the proband  called and no-call | For full context |
| **Father** | All 66 loci for the father | Cross-reference |
| **Mother** | All 66 loci for the mother | Cross-reference |
| **Legend** | Colour key | Share with clinical colleagues |

### Row Colour Coding

| Colour | Meaning |
|--------|---------|
| 🟣 Purple | De novo candidate  allele outside normal range, absent in both parents |
| 🟠 Orange | Population outlier HIGH  expansion above normal range |
| 🔵 Blue | Population outlier LOW  contraction below normal range |
| ⬜ Grey | Ambiguous call  nested/complex locus, treat with caution |
| 🟢 Green | PASS  within established normal repeat range |
| 🟡 Yellow | Low depth  called but coverage is low |
| ░ Light grey | No-call  insufficient coverage to genotype (REPCN = ./.) |

### De Novo vs Pathogenicity

> **Important:** `Allele_DeNovo = Yes` means the repeat was not observed in either parent  it is an **inheritance classification**, not a pathogenicity call.
>
> Clinical significance is determined by comparing the repeat count against the `PathogenicCutoff` column. A de novo repeat below the pathogenic threshold falls in the intermediate or reduced penetrance range and requires clinical interpretation in full context.
>
> Example: ATXN1 in LH0014A has 35 CAG repeats (de novo, above normal max of 32). The pathogenic cutoff for SCA1 is ≥39. This warrants clinical attention but is not a confirmed pathogenic expansion.

---

## REViewer  On-Demand Visualisation

[REViewer](https://github.com/Illumina/REViewer) generates read-level visualisations of STR loci using the realigned BAM produced by EH seeking mode. It is intentionally **not** integrated into the pipeline  use it on demand for specific candidate findings.

```bash
# Example: ATXN1 de novo candidate in LH0014A
REViewer \
    --reads     results/expansionhunter/LH0014A-RPT_S5/LH0014A-RPT_S5_realigned.sorted.bam \
    --vcf       results/expansionhunter/LH0014A-RPT_S5/LH0014A-RPT_S5.vcf \
    --reference /path/to/hg38.fasta \
    --catalog   /path/to/variant_catalog_v2.json \
    --locus     ATXN1 \
    --output-prefix results/reviewer/LH0014A-RPT_S5_ATXN1
```

Outputs: `*.svg` (open in browser), `*_metrics.tsv` (allele depth), `*_phasing.tsv` (diplotype string).

---

## Resource Configuration

The pipeline is configured for the RCSI HPC with a 30-CPU cap:

| Resource | Value | Notes |
|----------|-------|-------|
| Concurrent EH jobs | 25 | `params.slurm_queue_size` |
| CPUs per EH job | 1 | EH v5 is single-threaded |
| Peak total CPUs | 26 | 25 workers + 1 head job ≤ 30 cap |
| Memory per EH job | 4 GB | Retry bumps to 8 GB automatically |
| Aggregate step memory | 8 GB | openpyxl Excel generation |
| Sort temp directory | `$TMPDIR` | Node-local SSD  avoids NFS saturation |

Estimated runtime: 252 samples ÷ 25 concurrent × ~2h/sample ≈ 21h.

---

## Limitations

1. **WES coverage at STR loci is low**  ~40-50 of 65 loci are typically no-call per sample. WES capture kits do not enrich intronic/intergenic repeat loci. WGS would recover the missing loci at higher cost.

2. **Repeat count ≠ pathogenicity** : A HIGH outlier or de novo flag does not diagnose disease. Repeat counts must be compared against published pathogenic thresholds and interpreted in full clinical context, accounting for penetrance and expressivity.

3. **De novo requires a complete trio** : `DeNovo = Possible` (DUO) or blank (SINGLETON) should prompt parental sequencing before clinical interpretation.

4. **STRipy population reference**  Z-scores compare against 1KG Dragen dataset used by STRipy. Pathogenic cutoffs are from published literature and may be updated as evidence accumulates.

6. **Nested/complex loci (AmbiguousNestedCall)**  Some loci (e.g. TBX1) have complex architecture where EH detects a variant but cannot count the pathogenic motif. These are flagged grey in the Excel report and should not be interpreted without orthogonal validation.

---

## Repository Structure

```
.
├── main.nf                          # Main workflow
├── nextflow.config                  # Pipeline parameters and profiles
├── run_expansionhunter.slurm        # SLURM launcher
├── conf/
│   ├── base.config                  # Resource labels and check_max()
│   ├── slurm.config                 # SLURM executor settings
│   ├── local.config                 # Local execution (testing)
│   └── test.config                  # Minimal config for dry-run
├── modules/
│   ├── samtools.nf                  # SAMTOOLS_INDEX, EXTRACT_REGIONS, INDEX_REALIGNED
│   ├── expansionhunter.nf           # EXPANSIONHUNTER
│   └── manifest.nf                  # collectFile completion manifest
├── modules/local/stripy/
│   ├── fetch_locus_ref.nf           # STRIPY_FETCH_LOCUS_REF
│   ├── annotate_vcf.nf              # STRIPY_ANNOTATE_VCF
│   ├── build_sample_report.nf       # STRIPY_BUILD_SAMPLE_REPORT
│   └── aggregate_cohort_report.nf   # STRIPY_AGGREGATE_COHORT_REPORT
├── subworkflows/local/
│   └── stripy_postprocess.nf        # STRIPY_POSTPROCESS orchestration
└── scripts/
    ├── build_sample_report.py       # 48-column TSV generator
    ├── aggregate_cohort_report.py   # Per-family Excel generator
    ├── fetch_stripy_locus_ref.py    # STRipy /locus bootstrap
    ├── generate_regions_bed.py      # Padded repeat loci BED
 
```

---

## Citation

If you use this pipeline in your research, please cite:

- **ExpansionHunter:** Dolzhenko et al., (2019) https://doi.org/10.1093/bioinformatics/btz431
- **STRipy:** Halman A. et al., (2022) https://doi.org/10.1002/humu.24382
- **Nextflow:** Di Tommaso et al. (2017) https://doi.org/10.1038/nbt.3820

---

## Contact

**Deepak Bharti**  
Clinical Bioinformatician, RCSI  
deepakbharti@rcsi.ie

---

*Pipeline developed for the Epilepsy Project, RCSI  April 2026*

---
---
Please visit https://github.com/Deep2106/BUILD_SINGULARITY_IMAGES/tree/main for more information about containers used.
---
---
