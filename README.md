# Single-Cell eQTL Analysis of Glioma GWAS Risk Loci

Locus-level companion to the cell-stratified Mendelian Randomization (csMR)
analysis.  For each of 12 glioma GWAS lead SNPs, identifies which brain cell
type's regulatory architecture is most affected using single-nucleus eQTL data.

## Quick start

```bash
# 1. Install dependencies (into your existing Python environment)
cd /path/to/pipeline
pip install -r requirements.txt

# 2. Run everything
bash run_pipeline.sh

# All outputs land in ./output/
```

## Pipeline overview

```
 Step 0                Step 1                 Step 2
 00_download       01_expression          02_eqtl_lookup
 _bryois.sh        _profiling.py          .py
   │ (parallel         │ (CELLxGENE            │ (parallel across
   │  xargs -P)        │  Census query)        │  cell type × chrom,
   │                   │                       │  BH correction,
   │                   │                       │  permutation test)
   ▼                   ▼                       ▼
 data/bryois/     celltype_            eqtl_results_all.csv
 *.gz              expression.csv      eqtl_corrected.csv
                                       locus_summary.csv
                                       permutation_results.csv
   ┌───────────────────────────────────────┘
   │
   ▼  Step 3                Step 4               Step 5
 03_colocalization       04_visualization      05_compile_results
 .py (optional)          .py                   .py
   │ (parallel per       │                      │
   │  locus, needs       │                      │
   │  GWAS sumstats)     │                      │
   ▼                     ▼                      ▼
 coloc_results.csv    fig1–fig4.png         comprehensive_results.csv
 coloc_summary.csv                          interpretation.md
```

## Step-by-step

### Step 0 — Download Bryois eQTL data

```bash
bash 00_download_bryois.sh [N_JOBS]
```

Downloads Bryois et al. 2022 cell-type eQTL summary statistics from Zenodo
(~600 MB).  Parallelises with `xargs -P`.  Caches files — safe to re-run.

**Run this on a node with internet access** (e.g., login node).

### Step 1 — Expression profiling

```bash
python3 01_expression_profiling.py
```

Queries CZ CELLxGENE Census for per-cell-type expression of all 14 GWAS genes
across 8 brain cell types from normal human brain tissue.

**Requires internet** (S3 access).  ~8 GB RAM, 5–15 minutes.

#### Known issue: CELLxGENE Census query

Both `01_expression_profiling.py` and `01_laptop_query.py` have practical
problems with the current Census version (2025-11-08):

1. **`tiledbsoma` requires AVX2** CPU instructions and direct S3 access —
   often unavailable on cluster login/compute nodes.
2. **The original dataset ID** (`6f7fd0f1-...`, Allen Brain Cell Atlas DLPFC)
   **no longer exists** in Census 2025-11-08.
3. **Falling back to `tissue_general == 'brain'`** matches the entire Census
   brain corpus (tens of millions of cells), making the query impractically
   slow (~2+ hours, >14 GB RAM) even on a local machine.

**Recommended: use the pre-generated HPA expression CSV** provided with this
pipeline (`celltype_expression.csv`).  It contains mean expression (nCPM) from
Human Protein Atlas v24 single-nuclei brain RNA-seq for all 14 genes × 8 cell
types.  Place it in `output/` and run with `--skip-census`.

A future fix would narrow the Census query to a specific tissue (e.g.,
`tissue == 'cerebral cortex'`) or identify the correct replacement dataset ID.

### Step 2 — eQTL lookup (parallelised)

```bash
python3 02_eqtl_lookup.py [--cpus N]
```

For each GWAS locus, extracts Bryois eQTL results across all cell types.
Uses `ProcessPoolExecutor` across cell-type × chromosome files.

Key improvements over the original single-threaded script:

| Issue | Fix |
|---|---|
| No multiple testing correction | BH FDR within each locus across 6 cell types |
| Cherry-picks best of all proxies | Uses only the highest-r² proxy per locus |
| Pseudobulk competes in ranking | Excluded from "best cell type" (reference only) |
| Line-by-line scanning per SNP | Pre-builds lookup set, single-pass per file |
| No enrichment test | Permutation test (10,000 shuffles) |

**Runs offline.  Auto-detects CPUs.**

### Step 3 — Colocalization (optional)

```bash
python3 03_colocalization.py [--cpus N]
```

Runs Bayesian colocalization (coloc.abf, Giambartolomei et al. 2014) for each
locus × cell type pair.  Pure Python implementation — no R dependency.

**Requires GWAS summary statistics.**  Place them in `data/gwas/`:

```
data/gwas/gwas_sumstats.tsv.gz     # single file, all chromosomes
# OR
data/gwas/gwas_chr1.tsv.gz         # per-chromosome files
data/gwas/gwas_chr2.tsv.gz
...
```

Expected columns (tab-separated):

| Column | Description |
|--------|-------------|
| `snp_id` | rsID (e.g., rs12345) |
| `chr` | Chromosome (1–22, no "chr" prefix) |
| `pos` | Position (hg38) |
| `beta` | Effect size (log-OR for case-control) |
| `se` | Standard error of beta |
| `pvalue` | Association p-value |
| `eaf` | Effect allele frequency |

If no GWAS files are found, the script prints instructions and exits cleanly.
The rest of the pipeline still works.

### Step 4 — Visualization

```bash
python3 04_visualization.py
```

Generates publication-ready figures:

| Figure | Description |
|--------|-------------|
| `fig1_expression_dotplot.png` | 12 genes × 8 cell types (split by subtype) |
| `fig2_eqtl_heatmap.png` | Signed -log10(p) heatmap with BH-corrected stars |
| `fig3_locus_summary.png` | Best cell type per locus; hatching = fails BH |
| `fig4_coloc_heatmap.png` | PP.H4 heatmap (only if coloc ran) |

### Step 5 — Compile results

```bash
python3 05_compile_results.py
```

Merges all upstream outputs into `comprehensive_results.csv` (one row per locus)
and writes `interpretation.md` with project context and caveats.

## Master runner

```bash
bash run_pipeline.sh [OPTIONS]
```

| Option | Effect |
|--------|--------|
| `--cpus N` | Force N worker processes (default: auto-detect via `nproc`) |
| `--skip-download` | Skip step 0 (data already cached) |
| `--skip-census` | Skip step 1 (no internet, or expression CSV already exists) |
| `--skip-coloc` | Skip step 3 (no GWAS summary stats) |

### Typical cluster workflow

```bash
# On login node (has internet):
cd /path/to/pipeline
bash 00_download_bryois.sh

# Option A: run expression profiling on login node (if tiledbsoma works)
python3 01_expression_profiling.py

# Option B: run on laptop instead (if login node lacks AVX2 or S3 access)
#   On laptop:  python3 01_laptop_query.py
#   Then:       scp output/celltype_expression.csv you@cluster:.../pipeline/output/

# Submit compute job:
sbatch --cpus-per-task=16 --mem=16G --time=1:00:00 \
  /path/to/pipeline/run_pipeline.sh --skip-download --skip-census --cpus 16
```

## Resource requirements

| Resource | Step 0 | Step 1 | Steps 2–5 |
|----------|--------|--------|-----------|
| Network  | Zenodo | S3 (us-west-2) | None |
| RAM      | <1 GB  | ~8 GB  | ~4 GB |
| Disk     | ~600 MB | ~200 MB | ~50 MB |
| CPUs     | scales to 12 | 1 | scales to all |
| Runtime  | 3–10 min | 5–15 min | 2–5 min |

## Data sources

| Dataset | URL | Use |
|---------|-----|-----|
| CELLxGENE Census 2025-11-08 | https://chanzuckerberg.github.io/cellxgene-census/ | Per-cell-type expression |
| Bryois eQTL (Zenodo 7276971) | https://zenodo.org/records/7276971 | 8 cell types × 22 chromosomes |
| Bryois SNP positions | (same Zenodo record) | rsID → hg38 mapping |

## Configuration

All tuneable parameters live in `config.py`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CENSUS_VERSION` | `"2025-11-08"` | CELLxGENE snapshot (pin for reproducibility) |
| `BH_FDR` | 0.05 | Within-locus BH correction threshold |
| `COLOC_PP4` | 0.80 | PP.H4 threshold for colocalization |
| `N_PERM` | 10,000 | Permutations for enrichment test |

### Adding GWAS effect sizes for risk-allele direction

In `config.py`, populate the `gwas_beta`, `gwas_se`, and `gwas_eaf` fields
for each locus in `GWAS_LOCI`.  This enables the risk-allele concordance
column in the summary table.

```python
"PHLDB1": {
    "rsid": "rs12803321", "chr": 11, ...
    "gwas_beta": 0.25,          # log-OR from GWAS
    "gwas_se":   0.05,          # standard error
    "gwas_eaf":  0.42,          # effect allele frequency
},
```

## Known caveats

1. **rs55705857 (CCDC26)** is not genotyped in Bryois.  The highest-r²
   EUR proxy (rs143893586, r²=0.636) is used.  Effect sizes are attenuated.

2. **rs634537 (CDKN2A)** is not in Bryois.  rs613312 is used; the LD r²
   has not been formally calculated.  **Resolve before publication.**

3. **N=12 loci** is too small for formal group-level tests.  Fisher exact
   and permutation results are descriptive.

4. **Bryois sample size (N=196)** limits power for less-abundant cell types
   (OPC, microglia).

5. **Within-locus BH** is conservative for only 6 tests per locus —
   some real effects may not survive correction.

## Directory structure

```
pipeline/
├── config.py                  # Shared configuration (single source of truth)
├── 00_download_bryois.sh      # Parallel data download
├── 01_expression_profiling.py # CELLxGENE Census query (cluster)
├── 01_laptop_query.py         # CELLxGENE Census query (laptop workaround)
├── 02_eqtl_lookup.py          # Parallel eQTL extraction + BH correction
├── 03_colocalization.py       # coloc.abf (optional, needs GWAS data)
├── 04_visualization.py        # Publication-ready figures
├── 05_compile_results.py      # Master table + interpretation
├── run_pipeline.sh            # Orchestrator with CLI options
├── requirements.txt           # Python dependencies
├── README.md                  # This file
├── data/
│   ├── bryois/                # Downloaded eQTL files (auto-populated)
│   └── gwas/                  # User-supplied GWAS summary stats
└── output/                    # All results, figures, logs
```

## Reference

Bryois J et al. Cell-type-specific cis-eQTLs in eight human brain cell types
identify novel risk genes for psychiatric and neurological disorders.
*Nat Neurosci* 25, 1104–1112 (2022).
DOI: [10.1038/s41593-022-01128-z](https://doi.org/10.1038/s41593-022-01128-z)
