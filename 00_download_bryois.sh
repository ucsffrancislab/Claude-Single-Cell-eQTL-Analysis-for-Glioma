#!/usr/bin/env bash
set -euo pipefail
#
# 00_download_bryois.sh — Download Bryois et al. 2022 eQTL data from Zenodo.
#
# Parallelises file downloads with xargs -P, using all available CPUs.
# Run this BEFORE the analysis scripts.  Works on login nodes with internet.
#
# Usage:
#   cd /path/to/pipeline
#   bash 00_download_bryois.sh            # auto-detect CPUs
#   bash 00_download_bryois.sh 8          # force 8 parallel jobs
#

# ── Locate pipeline root ─────────────────────────────────────────────────────
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" \
        | awk '/Command=/{sub(/.*Command=/, ""); print $1}')")
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
DATA_DIR="${SCRIPT_DIR}/data/bryois"
mkdir -p "${DATA_DIR}"

# ---------- CPU detection ----------
if [[ $# -ge 1 ]]; then
    NJOBS="$1"
else
    NJOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
fi
# cap at 12 — Zenodo rate-limits aggressive downloaders
NJOBS=$(( NJOBS > 12 ? 12 : NJOBS ))
echo "=== Bryois eQTL data download (${NJOBS} parallel jobs) ==="
echo "    Target: ${DATA_DIR}"
echo ""

# ---------- Zenodo record metadata ----------
RECORD_ID="7276971"
ZENODO_API="https://zenodo.org/api/records/${RECORD_ID}"
LISTING="${DATA_DIR}/zenodo_files.json"

if [[ ! -f "${LISTING}" ]]; then
    echo "Fetching Zenodo record metadata..."
    curl -sL "${ZENODO_API}" -o "${LISTING}"
fi

# ---------- Build download list ----------
# Cell type file prefixes (must match config.py BRYOIS_ALL_PREFIX)
CT_PREFIXES=(
    "OPCs...COPs"
    "Astrocytes"
    "Excitatory.neurons"
    "Inhibitory.neurons"
    "Microglia"
    "Oligodendrocytes"
    "pb"
)

# Chromosomes containing our GWAS loci
CHROMS=(1 2 3 5 7 8 9 11 15 20)

DOWNLOAD_LIST="${DATA_DIR}/.download_urls.txt"
> "${DOWNLOAD_LIST}"

# SNP position file
SNP_FILE="snp_pos.txt.gz"
SNP_PATH="${DATA_DIR}/${SNP_FILE}"
if [[ ! -f "${SNP_PATH}" ]] || [[ $(stat -c%s "${SNP_PATH}" 2>/dev/null || stat -f%z "${SNP_PATH}" 2>/dev/null) -lt 1000 ]]; then
    SNP_URL=$(python3 -c "
import json
with open('${LISTING}') as f:
    rec = json.load(f)
files = {f['key']: f['links']['self'] for f in rec['files']}
print(files['${SNP_FILE}'])
")
    echo "${SNP_URL} ${SNP_PATH}" >> "${DOWNLOAD_LIST}"
fi

# eQTL files: cell_type x chromosome
for prefix in "${CT_PREFIXES[@]}"; do
    for chr in "${CHROMS[@]}"; do
        fname="${prefix}.${chr}.gz"
        fpath="${DATA_DIR}/${fname}"
        # Skip if already downloaded and > 1KB
        if [[ -f "${fpath}" ]] && [[ $(stat -c%s "${fpath}" 2>/dev/null || stat -f%z "${fpath}" 2>/dev/null) -gt 1000 ]]; then
            continue
        fi
        # Extract URL from Zenodo listing
        url=$(python3 -c "
import json, sys
with open('${LISTING}') as f:
    rec = json.load(f)
files = {f['key']: f['links']['self'] for f in rec['files']}
key = '${fname}'
if key in files:
    print(files[key])
else:
    print('MISSING', file=sys.stderr)
    sys.exit(1)
" 2>/dev/null) || { echo "  WARN: ${fname} not on Zenodo (skipped)"; continue; }
        echo "${url} ${fpath}" >> "${DOWNLOAD_LIST}"
    done
done

N_FILES=$(wc -l < "${DOWNLOAD_LIST}" | tr -d ' ')
if [[ "${N_FILES}" -eq 0 ]]; then
    echo "All files already cached — nothing to download."
    echo "=== Done ==="
    exit 0
fi

echo "Downloading ${N_FILES} files (${NJOBS} in parallel)..."
echo ""

# ---------- Parallel download ----------
# Each line: "URL DEST_PATH"
# curl -sL -o <dest> <url>  with retry
download_one() {
    local url="$1"
    local dest="$2"
    local fname
    fname=$(basename "${dest}")
    for attempt in 1 2 3; do
        if curl -sL --retry 2 --connect-timeout 30 -o "${dest}" "${url}"; then
            local size
            size=$(stat -c%s "${dest}" 2>/dev/null || stat -f%z "${dest}" 2>/dev/null)
            printf "  %-40s %6.1f MB\n" "${fname}" "$(echo "${size}/1048576" | bc -l 2>/dev/null || python3 -c "print(${size}/1048576:.1f)")"
            return 0
        fi
        echo "  RETRY ${attempt}/3: ${fname}"
        sleep $((attempt * 2))
    done
    echo "  FAILED: ${fname}"
    return 1
}
export -f download_one

cat "${DOWNLOAD_LIST}" | xargs -n 2 -P "${NJOBS}" bash -c 'download_one "$@"' _

# ---------- Verify ----------
echo ""
echo "=== Verification ==="
TOTAL_SIZE=0
N_PRESENT=0
N_EXPECTED=$(( ${#CT_PREFIXES[@]} * ${#CHROMS[@]} + 1 ))  # +1 for snp_pos
for prefix in "${CT_PREFIXES[@]}"; do
    for chr in "${CHROMS[@]}"; do
        fpath="${DATA_DIR}/${prefix}.${chr}.gz"
        if [[ -f "${fpath}" ]] && [[ $(stat -c%s "${fpath}" 2>/dev/null || stat -f%z "${fpath}" 2>/dev/null) -gt 1000 ]]; then
            fsize=$(stat -c%s "${fpath}" 2>/dev/null || stat -f%z "${fpath}" 2>/dev/null)
            TOTAL_SIZE=$((TOTAL_SIZE + fsize))
            N_PRESENT=$((N_PRESENT + 1))
        fi
    done
done
if [[ -f "${SNP_PATH}" ]]; then
    N_PRESENT=$((N_PRESENT + 1))
    fsize=$(stat -c%s "${SNP_PATH}" 2>/dev/null || stat -f%z "${SNP_PATH}" 2>/dev/null)
    TOTAL_SIZE=$((TOTAL_SIZE + fsize))
fi

echo "  Files present: ${N_PRESENT} / ${N_EXPECTED}"
echo "  Total size:    $(python3 -c "print(f'{${TOTAL_SIZE}/1e6:.0f}')") MB"

if [[ "${N_PRESENT}" -lt "${N_EXPECTED}" ]]; then
    echo ""
    echo "  WARNING: Some files missing. Re-run this script to retry."
    exit 1
fi

echo ""
echo "=== Download complete ==="
