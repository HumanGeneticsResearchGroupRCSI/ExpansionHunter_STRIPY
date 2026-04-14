#!/bin/bash
#SBATCH --job-name=nf_expansionhunter
#SBATCH --output=nf_expansionhunter_%j.out
#SBATCH --error=nf_expansionhunter_%j.err
#SBATCH --partition=<PART>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=deepakbharti@rcsi.ie

# ============================================================================
# EXPANSIONHUNTER NEXTFLOW PIPELINE LAUNCHER
# ============================================================================
# This is the HEAD/MASTER job that orchestrates all worker SLURM jobs.
# ExpansionHunter does NOT run inside this job.
#
# CPU BUDGET (30-CPU cap):
#   - This head job: 1 CPU
#   - Per EH worker: 1 CPU (single-threaded)
#   - Concurrent jobs: 25 (queueSize in slurm.config)
#   - Peak total: 26 CPUs ≤ 30 cap ✔
#
# MEMORY BUDGET:
#   - This head job: 4 GB
#   - Per EH worker: 4 GB (WES streaming mode)
#   - Peak total: ~104 GB (25 × 4 + 4)
#
# WALLTIME:
#   - 252 samples ÷ 25 concurrent × ~2h/sample ≈ 21h worst case
#   - Set to 48h for safety margin
# ============================================================================


# ============================================================================
# USER CONFIGURATION - Edit these paths
# ============================================================================

PROJECT_DIR="/path/to/TT"
PIPELINE_DIR="${PROJECT_DIR}/wes-expansionhunter-march2026"

PED="/path/to/TEST_DVG/results_v4/ped/master.ped"
BAM_DIR="/path/to/TEST_DVG/results_v4/bam/markdup"
REFERENCE="/path/toreference/hg38.fasta"
VARIANT_CATALOG="/path/toreference/variant_catalog_v2_extended_offtarget_added_chr.json"

OUTDIR="/path/to/TEST_DVG/results_v4/repeat_analysis"
BAM_SUFFIX=".md.bam"              # or '.markdup.recal.bam'
ANALYSIS_MODE="seeking"          # 'seeking' (hybrid mini-BAM, recommended) | 'streaming' (no realigned BAM)
REGIONS_BED="/path/toreference/repeat_regions.bed"
                                   # Generate once: python3 scripts/generate_regions_bed.py \
                                   #   --catalog ${VARIANT_CATALOG} --output ${REGIONS_BED}

# ── STRipy post-processing ────────────────────────────────────────────────
LOCUS_REF="/path/toreference/stripy_locus_reference.json"
HPO_FILE="/path/toresources/hpo_genes_to_phenotype.txt"
# Override per project: set to hpo_renal_terms.txt, hpo_cf_terms.txt etc.
# HPO_FILTER_FILE="/path/toresources/hpo_epilepsy_terms.txt"
COHORT_ID="LH_WES"

# ============================================================================
# LOAD MODULES
# ============================================================================


module purge
module load Java/17
# ============================================================================
# NEXTFLOW ENVIRONMENT (matching your working pattern)
# ============================================================================

# Work directory (NXF_WORK as env var, not -w flag)
export NXF_WORK="${PROJECT_DIR}/work2"
mkdir -p "${NXF_WORK}"

# Singularity image cache
export NXF_SINGULARITY_CACHEDIR="${PROJECT_DIR}/.singularity_cache"
mkdir -p "${NXF_SINGULARITY_CACHEDIR}"

# Java temp files → node-local SSD
# RCSI sets TMPDIR to /local/scratch automatically
export _JAVA_OPTIONS="-Djava.io.tmpdir=${TMPDIR:-/tmp}"

# JVM heap for Nextflow master
export NXF_OPTS="-Xms512m -Xmx3g"

# Cleaner log output
export NXF_ANSI_LOG=false

# ============================================================================
# SETUP DIRECTORIES
# ============================================================================

mkdir -p "${OUTDIR}/pipeline_info"


# ============================================================================
# PRE-FLIGHT CHECKS
# ============================================================================

echo "================================================================"
echo "  ExpansionHunter Pipeline - Pre-flight Checks"
echo "================================================================"
echo "  SLURM job ID    : ${SLURM_JOB_ID}"
echo "  Node            : $(hostname)"
echo "  Started         : $(date)"
echo ""
echo "  Pipeline dir    : ${PIPELINE_DIR}"
echo "  NXF_WORK        : ${NXF_WORK}"
echo "  NXF_SINGULARITY : ${NXF_SINGULARITY_CACHEDIR}"
echo "  TMPDIR          : ${TMPDIR:-/tmp (not set)}"
echo "  Output dir      : ${OUTDIR}"
echo "================================================================"

fail=0
for f in "${PED}" "${REFERENCE}" "${VARIANT_CATALOG}" "${BAM_DIR}" "${PIPELINE_DIR}/main.nf"        "${LOCUS_REF}" "${HPO_FILE}" "${REGIONS_BED}"; do
    if [[ ! -e "${f}" ]]; then
        echo "[ERROR] Not found: ${f}"
        fail=1
    else
        echo "[OK]    ${f}"
    fi
done

[[ ${fail} -eq 1 ]] && {
    echo ""
    echo "[ABORT] Fix missing files above"
    exit 1
}

SAMPLE_COUNT=$(grep -cE '^[^#[:space:]]' "${PED}" || true)

echo ""
echo "[INFO] Samples in PED     : ${SAMPLE_COUNT}"
echo "[INFO] BAM suffix         : ${BAM_SUFFIX}"
echo "[INFO] Analysis mode      : ${ANALYSIS_MODE}"
echo "[INFO] Concurrent jobs    : 25"
echo "[INFO] Est. runtime       : ~$(( (SAMPLE_COUNT / 25 + 1) * 2 ))h"
echo ""
echo "[INFO] Regions BED        : ${REGIONS_BED}"
echo "[INFO] STRipy locus ref   : ${LOCUS_REF}"
echo "[INFO] HPO file           : ${HPO_FILE}"
echo "[INFO] Cohort ID          : ${COHORT_ID}"
echo ""

# ============================================================================
# LAUNCH NEXTFLOW
# ============================================================================

echo "[$(date +%H:%M:%S)] Launching Nextflow pipeline..."
echo ""

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile          slurm \
    -resume \
    --ped             "${PED}" \
    --bam_dir         "${BAM_DIR}" \
    --fasta           "${REFERENCE}" \
    --variant_catalog "${VARIANT_CATALOG}" \
    --bam_suffix      "${BAM_SUFFIX}" \
    --analysis_mode   "${ANALYSIS_MODE}" \
    --outdir          "${OUTDIR}" \
    --stripy_locus_ref  "${LOCUS_REF}" \
    --hpo_file          "${HPO_FILE}" \
    --cohort_id         "${COHORT_ID}" \
    --regions_bed       "${REGIONS_BED}"

NXF_EXIT=$?

# ============================================================================
# POST-RUN SUMMARY
# ============================================================================

echo ""
echo "================================================================"
if [[ ${NXF_EXIT} -eq 0 ]]; then
    echo "  STATUS   : SUCCESS ✔"
    echo ""
    echo "  Results  : ${OUTDIR}"
    echo ""
    echo "  Clean up work cache:"
    echo "    nextflow clean -f"
    echo "    rm -rf ${NXF_WORK}"
else
    echo "  STATUS   : FAILED ✘ (exit: ${NXF_EXIT})"
    echo ""
    echo "  Diagnose:"
    echo "    grep -E 'ERROR|FAILED' logs/nf_expansionhunter_${SLURM_JOB_ID}.out"
    echo "    nextflow log last -f name,status,exit,workdir"
    echo ""
    echo "  Re-run:"
    echo "    sbatch run_expansionhunter.slurm"
fi
echo "  Finished : $(date)"
echo "================================================================"

exit ${NXF_EXIT}
