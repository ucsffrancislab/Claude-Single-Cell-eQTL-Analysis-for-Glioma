#!/usr/bin/env bash
set -euo pipefail
#
# run_pipeline.sh — Master orchestrator for the glioma sc-eQTL pipeline.
#
# Usage:
#   cd /path/to/pipeline
#   bash run_pipeline.sh [--cpus N] [--skip-download] [--skip-census] [--skip-coloc]
#
# The pipeline runs 6 steps in order:
#   Step 0: Download Bryois data (parallel, needs internet)
#   Step 1: CELLxGENE expression profiling (needs internet / S3)
#   Step 2: eQTL lookup (CPU-parallel, offline OK)
#   Step 3: Colocalization (optional, needs GWAS sumstats)
#   Step 4: Visualization (offline)
#   Step 5: Compile results & interpretation (offline)
#

# ── Locate pipeline root (where this script lives) ──────────────────────────
# Under SLURM, sbatch copies the script to a spool directory so dirname "$0"
# points to the wrong place.  Use scontrol to recover the original path.
# Outside SLURM (interactive/local), dirname "$0" works fine.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" \
        | awk '/Command=/{sub(/.*Command=/, ""); print $1}')")
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
export SCEQTL_PIPELINE_DIR="${SCRIPT_DIR}"

# Work directory: --work-dir, or SCEQTL_WORK_DIR, or current directory
if [[ -n "${WORK_DIR}" ]]; then
    export SCEQTL_WORK_DIR="$(cd "${WORK_DIR}" && pwd)"
elif [[ -z "${SCEQTL_WORK_DIR:-}" ]]; then
    export SCEQTL_WORK_DIR="$(pwd)"
fi
cd "${SCEQTL_WORK_DIR}"

# Ensure Python can find config.py and other pipeline modules in SCRIPT_DIR
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

# ---------- Parse arguments ----------
CPUS=0
SKIP_DOWNLOAD=false
SKIP_CENSUS=false
SKIP_COLOC=false

WORK_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpus)          CPUS="$2"; shift 2 ;;
        --work-dir)      WORK_DIR="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        --skip-census)   SKIP_CENSUS=true; shift ;;
        --skip-coloc)    SKIP_COLOC=true; shift ;;
        -h|--help)
            echo "Usage: bash run_pipeline.sh [--work-dir DIR] [--cpus N] [--skip-download] [--skip-census] [--skip-coloc]"
            echo ""
            echo "  --work-dir DIR   Data/output directory (default: current directory)"
            echo "  --cpus N         Worker processes (default: auto-detect)"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------- CPU detection ----------
if [[ "${CPUS}" -eq 0 ]]; then
    # SLURM-aware: use allocated CPUs, not total node CPUs
    CPUS=${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || 4)}
fi
echo "============================================================"
echo " Glioma Single-Cell eQTL Pipeline"
echo " $(date)"
echo " CPUs: ${CPUS}"
echo " Scripts:      ${SCRIPT_DIR}
 Work dir:     ${SCEQTL_WORK_DIR}"
echo "============================================================"
echo ""

# ---------- Dependency check ----------
echo "--- Checking dependencies ---"
MISSING=()

python3 -c "import pandas"      2>/dev/null || MISSING+=("pandas")
python3 -c "import numpy"       2>/dev/null || MISSING+=("numpy")
python3 -c "import scipy"       2>/dev/null || MISSING+=("scipy")
python3 -c "import matplotlib"  2>/dev/null || MISSING+=("matplotlib")
python3 -c "import statsmodels" 2>/dev/null || MISSING+=("statsmodels")
python3 -c "import requests"    2>/dev/null || MISSING+=("requests")

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "ERROR: Missing Python packages: ${MISSING[*]}"
    echo "Install with: pip install -r requirements.txt"
    exit 1
fi
echo "  Core packages OK"

# Check optional Census dependencies
if [[ "${SKIP_CENSUS}" == "false" ]]; then
    python3 -c "import cellxgene_census" 2>/dev/null || {
        echo "WARNING: cellxgene_census not installed — step 01 will fail."
        echo "  Install with: pip install cellxgene-census tiledbsoma"
        echo "  Or re-run with --skip-census to skip expression profiling."
        echo ""
    }
fi

echo ""

# ---------- Step 0: Download ----------
STEP_START=$(date +%s)
if [[ "${SKIP_DOWNLOAD}" == "true" ]]; then
    echo "=== Step 0: Download (SKIPPED) ==="
else
    echo "=== Step 0: Download Bryois data ==="
    bash "${SCRIPT_DIR}/00_download_bryois.sh" "${CPUS}"
fi
echo "  Step 0: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ---------- Step 1: Expression profiling ----------
STEP_START=$(date +%s)
if [[ "${SKIP_CENSUS}" == "true" ]]; then
    echo "=== Step 1: Expression profiling (SKIPPED) ==="
else
    echo "=== Step 1: Expression profiling (CELLxGENE Census) ==="
    python3 "${SCRIPT_DIR}/01_expression_profiling.py"
fi
echo "  Step 1: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ---------- Step 2: eQTL lookup ----------
STEP_START=$(date +%s)
echo "=== Step 2: Parallel eQTL lookup ==="
python3 "${SCRIPT_DIR}/02_eqtl_lookup.py" --cpus "${CPUS}"
echo "  Step 2: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ---------- Step 3: Colocalization ----------
STEP_START=$(date +%s)
if [[ "${SKIP_COLOC}" == "true" ]]; then
    echo "=== Step 3: Colocalization (SKIPPED) ==="
else
    echo "=== Step 3: Colocalization ==="
    python3 "${SCRIPT_DIR}/03_colocalization.py" --cpus "${CPUS}"
fi
echo "  Step 3: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ---------- Step 4: Visualization ----------
STEP_START=$(date +%s)
echo "=== Step 4: Visualization ==="
python3 "${SCRIPT_DIR}/04_visualization.py"
echo "  Step 4: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ---------- Step 5: Compile ----------
STEP_START=$(date +%s)
echo "=== Step 5: Compile results & interpretation ==="
python3 "${SCRIPT_DIR}/05_compile_results.py"
echo "  Step 5: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ---------- Summary ----------
echo "============================================================"
echo " Pipeline complete!"
echo " $(date)"
echo ""
echo " Outputs in: ${SCEQTL_WORK_DIR}/output/"
echo ""
ls -lh "${SCEQTL_WORK_DIR}"/output/*.csv "${SCEQTL_WORK_DIR}"/output/*.png "${SCEQTL_WORK_DIR}"/output/*.md "${SCEQTL_WORK_DIR}"/output/*.log 2>/dev/null || true
echo ""
echo " Start with: output/interpretation.md"
echo "============================================================"
