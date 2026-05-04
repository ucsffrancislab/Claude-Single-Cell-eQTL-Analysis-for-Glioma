#!/usr/bin/env python3
"""
03_colocalization.py
--------------------
Bayesian colocalization (coloc.abf) for each (locus x cell_type) pair.

Implements the Giambartolomei et al. 2014 approximate Bayes factor method
in pure Python — no R dependency.  Parallelised across loci with
multiprocessing.

REQUIRES
--------
  data/gwas/<subtype>_chr<N>.tsv.gz   (or a single all-chroms file)
      Columns: snp_id, chr, pos, beta, se, pvalue, eaf
      (log-OR scale for case-control; beta scale for quantitative)

  data/bryois/<prefix>.<chr>.gz       (downloaded by 00_download_bryois.sh)

If GWAS summary stats are not present, this script prints instructions and
exits cleanly — the rest of the pipeline can still run.

Output
------
  output/coloc_results.csv   — PP.H0–H4 per (locus, cell_type)
  output/coloc_summary.csv   — per-locus best PP.H4 and interpretation
"""

import os, sys, gzip, argparse, glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.stats import norm as sp_norm

from config import (
    GWAS_LOCI, GWAS_DIR, BRYOIS_DIR, OUTPUT_DIR,
    BRYOIS_CT_PREFIX, BRYOIS_N, NEEDED_CHROMS,
    RANKING_CELL_TYPES, COLOC_PP4,
    snps_for_locus, genes_for_locus, setup_logging,
)

log = setup_logging("03-coloc", OUTPUT_DIR / "03_colocalization.log")

# ── coloc.abf in Python ─────────────────────────────────────────────────────
#
# Reference: Giambartolomei et al. PLoS Genet 2014;10(5):e1004383
# Mirrors logic from R coloc::coloc.abf (v5.2+)

def _logsumexp(x):
    """Numerically stable log-sum-exp."""
    mx = np.max(x)
    return mx + np.log(np.sum(np.exp(x - mx)))


def wakefield_abf(beta, varbeta, prior_var):
    """Wakefield approximate Bayes factor (log scale).

    lABF = 0.5 * [ log(V/(V+W)) + z^2 * W/(V+W) ]
    V = varbeta, W = prior_var, z = beta / sqrt(V)
    """
    V = varbeta
    W = prior_var
    z2 = (beta ** 2) / V
    r = W / (V + W)
    return 0.5 * (np.log(1 - r) + z2 * r)


def coloc_abf(beta1, varbeta1, beta2, varbeta2,
              type1="cc", type2="quant",
              p1=1e-4, p2=1e-4, p12=1e-5,
              prior_var1=None, prior_var2=None):
    """
    Run coloc.abf for two traits sharing the same SNP set.

    Parameters
    ----------
    beta1, varbeta1 : array-like  (GWAS)
    beta2, varbeta2 : array-like  (eQTL)
    type1 : 'cc' (case-control) or 'quant'
    type2 : 'cc' or 'quant'
    p1, p2, p12 : per-SNP prior probabilities

    Returns
    -------
    dict with PP.H0 .. PP.H4, nsnps
    """
    beta1 = np.asarray(beta1, dtype=float)
    beta2 = np.asarray(beta2, dtype=float)
    varbeta1 = np.asarray(varbeta1, dtype=float)
    varbeta2 = np.asarray(varbeta2, dtype=float)

    assert len(beta1) == len(beta2), "SNP count mismatch"
    nsnps = len(beta1)

    # Default prior variances (coloc defaults)
    if prior_var1 is None:
        prior_var1 = 0.04 if type1 == "cc" else 0.15 ** 2
    if prior_var2 is None:
        prior_var2 = 0.04 if type2 == "cc" else 0.15 ** 2

    lABF1 = wakefield_abf(beta1, varbeta1, prior_var1)
    lABF2 = wakefield_abf(beta2, varbeta2, prior_var2)

    # Hypothesis likelihoods (log scale)
    lH0 = 0.0
    lH1 = np.log(p1)  + _logsumexp(lABF1)
    lH2 = np.log(p2)  + _logsumexp(lABF2)
    lH3 = np.log(p1)  + np.log(p2) + _logsumexp(lABF1) + _logsumexp(lABF2)
    lH4 = np.log(p12) + _logsumexp(lABF1 + lABF2)

    all_h = np.array([lH0, lH1, lH2, lH3, lH4])
    denom = _logsumexp(all_h)
    pp = np.exp(all_h - denom)

    return {
        "nsnps":  nsnps,
        "PP.H0":  round(pp[0], 4),
        "PP.H1":  round(pp[1], 4),
        "PP.H2":  round(pp[2], 4),
        "PP.H3":  round(pp[3], 4),
        "PP.H4":  round(pp[4], 4),
    }


# ── Data loading helpers ─────────────────────────────────────────────────────

def load_gwas_region(chrom, pos_center, window=500_000):
    """Load GWAS summary stats within ±window of pos_center.

    Looks for:  data/gwas/*chr{chrom}*.tsv.gz  OR  data/gwas/gwas_sumstats.tsv.gz
    Expected columns: snp_id, chr, pos, beta, se, pvalue, eaf
    """
    patterns = [
        GWAS_DIR / f"*chr{chrom}*.tsv.gz",
        GWAS_DIR / f"*chr{chrom}*.tsv",
        GWAS_DIR / "gwas_sumstats.tsv.gz",
        GWAS_DIR / "gwas_sumstats.tsv",
    ]
    for pat in patterns:
        hits = sorted(glob.glob(str(pat)))
        if hits:
            df = pd.read_csv(hits[0], sep="\t",
                             dtype={"chr": str, "snp_id": str})
            df["chr"] = df["chr"].astype(str).str.replace("chr", "")
            df = df[df["chr"] == str(chrom)]
            df = df[(df["pos"] >= pos_center - window)
                    & (df["pos"] <= pos_center + window)]
            return df
    return None


def load_bryois_region(ct_prefix, chrom, gene_symbol, snp_ids_in_region):
    """Load all Bryois eQTL rows for a gene on a chromosome,
    restricted to SNPs in snp_ids_in_region (set)."""
    fpath = BRYOIS_DIR / f"{ct_prefix}.{chrom}.gz"
    if not fpath.exists():
        return None
    rows = []
    with gzip.open(fpath, "rt") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 5:
                continue
            sym = parts[0].split("_")[0]
            if sym != gene_symbol:
                continue
            snp = parts[1]
            if snp not in snp_ids_in_region:
                continue
            rows.append({
                "snp_id":  snp,
                "beta":    float(parts[4]),
                "pvalue":  float(parts[3]),
            })
    if not rows:
        return None
    df = pd.DataFrame(rows)
    # Derive SE from beta and p-value:  z = beta/SE → SE = |beta|/|z|
    z = sp_norm.isf(df["pvalue"] / 2)        # two-sided
    z = np.where(z == 0, 1e-10, z)           # guard against p=1
    df["se"] = np.abs(df["beta"]) / np.abs(z)
    df.loc[df["se"] == 0, "se"] = 1e-10      # guard
    return df


def _get_snp_positions():
    """Load Bryois snp_pos.txt.gz → dict {snp_id: (chr, pos)}."""
    path = BRYOIS_DIR / "snp_pos.txt.gz"
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", compression="gzip",
                     dtype={"chr": str})
    return {r["snp"]: (r["chr"], r["pos"]) for _, r in df.iterrows()}


# ── Per-locus colocalization (runs in worker) ────────────────────────────────

def _coloc_one_locus(args):
    """Run coloc for one locus across all ranking cell types."""
    locus, snp_positions = args
    cfg = GWAS_LOCI[locus]
    chrom = cfg["chr"]
    gene = locus                    # primary gene
    lead_snp = cfg["rsid"]

    # Get genomic position of lead SNP
    pos_info = snp_positions.get(lead_snp)
    if pos_info is None:
        # Try proxies
        for px_id, _ in cfg.get("proxies", []):
            pos_info = snp_positions.get(px_id)
            if pos_info:
                break
    if pos_info is None:
        return [{"locus": locus, "_skip": "lead SNP position unknown"}]

    _, pos_center = pos_info
    pos_center = int(pos_center)

    # Load GWAS region
    gwas = load_gwas_region(chrom, pos_center)
    if gwas is None or gwas.empty:
        return [{"locus": locus, "_skip": "no GWAS data for region"}]

    gwas_snps = set(gwas["snp_id"].values)
    results = []
    for ct_name, ct_prefix in BRYOIS_CT_PREFIX.items():
        eqtl = load_bryois_region(ct_prefix, chrom, gene, gwas_snps)
        if eqtl is None or eqtl.empty:
            results.append({
                "locus": locus, "gene": gene, "cell_type": ct_name,
                "nsnps": 0, "PP.H0": None, "PP.H1": None,
                "PP.H2": None, "PP.H3": None, "PP.H4": None,
                "note": "no overlapping SNPs",
            })
            continue

        # Merge on SNP
        merged = eqtl.merge(gwas[["snp_id", "beta", "se"]],
                            on="snp_id", suffixes=("_eqtl", "_gwas"))
        if len(merged) < 10:
            results.append({
                "locus": locus, "gene": gene, "cell_type": ct_name,
                "nsnps": len(merged),
                "PP.H0": None, "PP.H1": None, "PP.H2": None,
                "PP.H3": None, "PP.H4": None,
                "note": f"<10 shared SNPs ({len(merged)})",
            })
            continue

        pp = coloc_abf(
            beta1=merged["beta_gwas"].values,
            varbeta1=(merged["se_gwas"].values ** 2),
            beta2=merged["beta_eqtl"].values,
            varbeta2=(merged["se_eqtl"].values ** 2),
            type1="cc", type2="quant",
        )
        pp["locus"] = locus
        pp["gene"] = gene
        pp["cell_type"] = ct_name
        pp["note"] = ""
        results.append(pp)

    return results


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Colocalization analysis")
    parser.add_argument("--cpus", type=int, default=0,
                        help="Worker processes (0 = auto-detect)")
    args = parser.parse_args()
    ncpus = args.cpus if args.cpus > 0 else (
        len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity')
        else os.cpu_count() or 4
    )

    # ── Check for GWAS data ──
    gwas_files = list(GWAS_DIR.glob("*.tsv*"))
    if not gwas_files:
        log.info("=" * 72)
        log.info("GWAS summary statistics not found — skipping colocalization")
        log.info("=" * 72)
        log.info("")
        log.info("To enable colocalization, place GWAS summary stats in:")
        log.info(f"  {GWAS_DIR}/")
        log.info("")
        log.info("Expected format (tab-separated, optionally gzipped):")
        log.info("  snp_id  chr  pos  beta  se  pvalue  eaf")
        log.info("")
        log.info("Accepted filenames:")
        log.info("  gwas_sumstats.tsv.gz    (single file, all chroms)")
        log.info("  gwas_chr1.tsv.gz        (per-chromosome files)")
        log.info("")
        log.info("The pipeline continues without colocalization.")

        # Write a stub output so downstream scripts don't break
        stub = pd.DataFrame(columns=["locus", "gene", "cell_type",
                                      "PP.H4", "note"])
        stub.to_csv(OUTPUT_DIR / "coloc_results.csv", index=False)
        return

    log.info("=" * 72)
    log.info("Colocalization analysis (coloc.abf)")
    log.info("=" * 72)
    log.info(f"  GWAS files: {[f.name for f in gwas_files]}")
    log.info(f"  Using {ncpus} workers")

    # Load SNP positions for region definition
    snp_positions = _get_snp_positions()
    log.info(f"  {len(snp_positions):,} SNP positions loaded")

    # Run per-locus in parallel
    tasks = [(locus, snp_positions) for locus in GWAS_LOCI]
    all_results = []
    with ProcessPoolExecutor(max_workers=ncpus) as pool:
        futures = {pool.submit(_coloc_one_locus, t): t[0] for t in tasks}
        for fut in as_completed(futures):
            locus = futures[fut]
            rows = fut.result()
            for r in rows:
                if "_skip" in r:
                    log.info(f"  [{locus}] skipped: {r['_skip']}")
                else:
                    all_results.append(r)

    coloc_df = pd.DataFrame(all_results)
    out = OUTPUT_DIR / "coloc_results.csv"
    coloc_df.to_csv(out, index=False)
    log.info(f"  Saved {out} ({len(coloc_df)} rows)")

    # Summary: best PP.H4 per locus
    if "PP.H4" in coloc_df.columns and coloc_df["PP.H4"].notna().any():
        valid = coloc_df.dropna(subset=["PP.H4"])
        idx = valid.groupby("locus")["PP.H4"].idxmax()
        summary = valid.loc[idx, ["locus", "gene", "cell_type",
                                   "nsnps", "PP.H4"]].copy()
        summary["colocalized"] = summary["PP.H4"] >= COLOC_PP4
        summary = summary.sort_values("PP.H4", ascending=False)
        sout = OUTPUT_DIR / "coloc_summary.csv"
        summary.to_csv(sout, index=False)
        log.info(f"  Saved {sout}")
        log.info(f"\n{summary.to_string(index=False)}\n")
    else:
        log.info("  No valid coloc results to summarise.")

    log.info("Done.")


if __name__ == "__main__":
    main()
