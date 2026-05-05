#!/usr/bin/env bash
set -euo pipefail
#
# run_pipeline.sh — Master orchestrator for the glioma sc-eQTL pipeline.
#
# Usage:
#   bash run_pipeline.sh [OPTIONS]
#
#   --work-dir DIR   Data/output directory (default: current directory)
#   --venv DIR       Activate Python venv at DIR before running
#   --cpus N         Worker processes (default: auto-detect from SLURM)
#   --skip-download  Skip Bryois data download
#   --skip-census    Skip CELLxGENE expression profiling
#   --skip-coloc     Skip colocalization
#

# ── Locate pipeline root (where this script lives) ──────────────────────────
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" \
        | awk '/Command=/{sub(/.*Command=/, ""); print $1}')")
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
export SCEQTL_PIPELINE_DIR="${SCRIPT_DIR}"

# ── Parse arguments ──────────────────────────────────────────────────────────
CPUS=0
SKIP_DOWNLOAD=false
SKIP_CENSUS=false
SKIP_COLOC=false
WORK_DIR=""
VENV=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpus)          CPUS="$2"; shift 2 ;;
        --work-dir)      WORK_DIR="$2"; shift 2 ;;
        --venv)          VENV="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        --skip-census)   SKIP_CENSUS=true; shift ;;
        --skip-coloc)    SKIP_COLOC=true; shift ;;
        -h|--help)
            echo "Usage: bash run_pipeline.sh [OPTIONS]"
            echo ""
            echo "  --work-dir DIR   Data/output directory (default: current directory)"
            echo "  --venv DIR       Activate Python venv at DIR before running"
            echo "  --cpus N         Worker processes (default: auto-detect from SLURM)"
            echo "  --skip-download  Skip Bryois data download"
            echo "  --skip-census    Skip CELLxGENE expression profiling"
            echo "  --skip-coloc     Skip colocalization"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Activate venv (if requested) ────────────────────────────────────────────
if [[ -n "${VENV}" ]]; then
    if [[ -f "${VENV}/bin/activate" ]]; then
        source "${VENV}/bin/activate"
    else
        echo "ERROR: venv not found at ${VENV}/bin/activate"
        exit 1
    fi
fi

# ── Work directory ───────────────────────────────────────────────────────────
if [[ -n "${WORK_DIR}" ]]; then
    export SCEQTL_WORK_DIR="$(cd "${WORK_DIR}" && pwd)"
elif [[ -z "${SCEQTL_WORK_DIR:-}" ]]; then
    export SCEQTL_WORK_DIR="$(pwd)"
fi
cd "${SCEQTL_WORK_DIR}"

# Ensure Python can find config.py and other pipeline modules
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

# ── CPU detection ────────────────────────────────────────────────────────────
if [[ "${CPUS}" -eq 0 ]]; then
    CPUS=${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || 4)}
fi

echo "============================================================"
echo " Glioma Single-Cell eQTL Pipeline"
echo " $(date)"
echo " CPUs:        ${CPUS}"
echo " Scripts:     ${SCRIPT_DIR}"
echo " Work dir:    ${SCEQTL_WORK_DIR}"
echo "============================================================"
echo ""

# ── Dependency check ─────────────────────────────────────────────────────────
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
    echo "Install with: pip install -r ${SCRIPT_DIR}/requirements.txt"
    exit 1
fi
echo "  Core packages OK"

if [[ "${SKIP_CENSUS}" == "false" ]]; then
    python3 -c "import cellxgene_census" 2>/dev/null || {
        echo "WARNING: cellxgene_census not installed — step 01 will fail."
        echo "  Install with: pip install cellxgene-census tiledbsoma"
        echo "  Or re-run with --skip-census"
        echo ""
    }
fi

echo ""

# ── Step 0: Download ─────────────────────────────────────────────────────────
STEP_START=$(date +%s)
if [[ "${SKIP_DOWNLOAD}" == "true" ]]; then
    echo "=== Step 0: Download (SKIPPED) ==="
else
    echo "=== Step 0: Download Bryois data ==="
    bash "${SCRIPT_DIR}/00_download_bryois.sh" "${CPUS}"
fi
echo "  Step 0: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ── Step 1: Expression profiling ─────────────────────────────────────────────
STEP_START=$(date +%s)
if [[ "${SKIP_CENSUS}" == "true" ]]; then
    echo "=== Step 1: Expression profiling (SKIPPED) ==="
else
    echo "=== Step 1: Expression profiling (CELLxGENE Census) ==="
    python3 "${SCRIPT_DIR}/01_expression_profiling.py"
fi
echo "  Step 1: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ── Step 2: eQTL lookup ─────────────────────────────────────────────────────
STEP_START=$(date +%s)
echo "=== Step 2: Parallel eQTL lookup ==="
python3 "${SCRIPT_DIR}/02_eqtl_lookup.py" --cpus "${CPUS}"
echo "  Step 2: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ── Step 3: Colocalization ───────────────────────────────────────────────────
STEP_START=$(date +%s)
if [[ "${SKIP_COLOC}" == "true" ]]; then
    echo "=== Step 3: Colocalization (SKIPPED) ==="
else
    echo "=== Step 3: Colocalization ==="
    python3 "${SCRIPT_DIR}/03_colocalization.py" --cpus "${CPUS}"
fi
echo "  Step 3: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ── Step 4: Visualization ───────────────────────────────────────────────────
STEP_START=$(date +%s)
echo "=== Step 4: Visualization ==="
python3 "${SCRIPT_DIR}/04_visualization.py"
echo "  Step 4: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ── Step 5: Compile ─────────────────────────────────────────────────────────
STEP_START=$(date +%s)
echo "=== Step 5: Compile results & interpretation ==="
python3 "${SCRIPT_DIR}/05_compile_results.py"
echo "  Step 5: $(( $(date +%s) - STEP_START ))s elapsed"
echo ""

# ── Summary ──────────────────────────────────────────────────────────────────
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
