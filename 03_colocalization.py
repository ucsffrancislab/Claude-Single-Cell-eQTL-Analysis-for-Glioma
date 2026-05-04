#!/usr/bin/env python3
"""
03_colocalization.py
--------------------
Bayesian colocalization (coloc.abf) for each (locus x cell_type) pair.

Implements the Giambartolomei et al. 2014 approximate Bayes factor method
in pure Python — no R dependency.  Parallelised across loci with
multiprocessing.

GWAS FORMAT
-----------
Reads your meta-analysis summary stats as-is (no reformatting needed).
Expected columns: CHR, BP, SNP, A1, A2, A1_FREQ, BETA, SE, P
(column names are mapped via GWAS_COLUMN_MAP in config.py)

Matching between GWAS and Bryois eQTL is by genomic position (chr + bp),
since GWAS uses chr:pos:a1:a2 IDs while Bryois uses rsIDs.

Subtype-matched: IDH-mut loci use the IDH-mut GWAS file, IDH-wt loci
use the IDH-wt GWAS file (configured in GWAS_SUBTYPE_FILES in config.py).

Output
------
  output/coloc_results.csv   — PP.H0–H4 per (locus, cell_type)
  output/coloc_summary.csv   — per-locus best PP.H4 and interpretation
"""

import os, sys, gzip, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.stats import norm as sp_norm

from config import (
    GWAS_LOCI, GWAS_DIR, BRYOIS_DIR, OUTPUT_DIR,
    BRYOIS_CT_PREFIX, BRYOIS_N, NEEDED_CHROMS,
    RANKING_CELL_TYPES, COLOC_PP4,
    GWAS_SUBTYPE_FILES, GWAS_COLUMN_MAP, GWAS_BUILD,
    snps_for_locus, genes_for_locus, setup_logging,
)

log = setup_logging("03-coloc", OUTPUT_DIR / "03_colocalization.log")


# ── coloc.abf in Python ─────────────────────────────────────────────────────

def _logsumexp(x):
    mx = np.max(x)
    return mx + np.log(np.sum(np.exp(x - mx)))


def wakefield_abf(beta, varbeta, prior_var):
    V = varbeta
    W = prior_var
    z2 = (beta ** 2) / V
    r = W / (V + W)
    return 0.5 * (np.log(1 - r) + z2 * r)


def coloc_abf(beta1, varbeta1, beta2, varbeta2,
              type1="cc", type2="quant",
              p1=1e-4, p2=1e-4, p12=1e-5,
              prior_var1=None, prior_var2=None):
    beta1 = np.asarray(beta1, dtype=float)
    beta2 = np.asarray(beta2, dtype=float)
    varbeta1 = np.asarray(varbeta1, dtype=float)
    varbeta2 = np.asarray(varbeta2, dtype=float)

    assert len(beta1) == len(beta2), "SNP count mismatch"
    nsnps = len(beta1)

    if prior_var1 is None:
        prior_var1 = 0.04 if type1 == "cc" else 0.15 ** 2
    if prior_var2 is None:
        prior_var2 = 0.04 if type2 == "cc" else 0.15 ** 2

    lABF1 = wakefield_abf(beta1, varbeta1, prior_var1)
    lABF2 = wakefield_abf(beta2, varbeta2, prior_var2)

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

def _get_snp_positions():
    """Load Bryois snp_pos.txt.gz -> dict {snp_id: (chr, pos)} and
    reverse dict {(chr, pos): snp_id} for position-based matching."""
    path = BRYOIS_DIR / "snp_pos.txt.gz"
    if not path.exists():
        return {}, {}
    df = pd.read_csv(path, sep="\t", compression="gzip",
                     dtype={"chr": str, "pos": int, "snp": str})
    rsid_to_pos = {}
    pos_to_rsid = {}
    for _, r in df.iterrows():
        chrom = str(r["chr"]).replace("chr", "")
        pos = int(r["pos"])
        rsid_to_pos[r["snp"]] = (chrom, pos)
        # Multiple rsIDs can map to same position; keep first
        key = (chrom, pos)
        if key not in pos_to_rsid:
            pos_to_rsid[key] = r["snp"]
    return rsid_to_pos, pos_to_rsid


def load_gwas_region(subtype, chrom, pos_center, window=500_000):
    """Load subtype-matched GWAS summary stats within ±window of pos_center.

    Reads the file in its original format and maps columns via GWAS_COLUMN_MAP.
    """
    gwas_filename = GWAS_SUBTYPE_FILES.get(subtype)
    if not gwas_filename:
        return None

    gwas_path = GWAS_DIR / gwas_filename
    if not gwas_path.exists():
        return None

    # Read only the columns we need (saves memory on 1.3 GB files)
    needed_cols = list(GWAS_COLUMN_MAP.keys())
    try:
        df = pd.read_csv(gwas_path, sep="\t", usecols=needed_cols,
                         dtype={"CHR": str, "BP": int})
    except Exception as e:
        log.error(f"  Failed to read {gwas_path}: {e}")
        return None

    # Rename to internal names
    df = df.rename(columns=GWAS_COLUMN_MAP)

    # Filter to region
    df["chr"] = df["chr"].astype(str).str.replace("chr", "")
    df = df[df["chr"] == str(chrom)]
    df = df[(df["pos"] >= pos_center - window) &
            (df["pos"] <= pos_center + window)]

    return df


def load_bryois_region(ct_prefix, chrom, gene_symbol, positions_in_region):
    """Load all Bryois eQTL rows for a gene on a chromosome,
    restricted to SNPs at positions in positions_in_region (set of (chr, pos)).

    Returns DataFrame with columns: snp_id, pos, beta, pvalue, se
    """
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
            rows.append({
                "snp_id":  snp,
                "beta":    float(parts[4]),
                "pvalue":  float(parts[3]),
            })

    if not rows:
        return None

    df = pd.DataFrame(rows)
    # Derive SE from beta and p-value
    z = sp_norm.isf(df["pvalue"] / 2)
    z = np.where(z == 0, 1e-10, z)
    df["se"] = np.abs(df["beta"]) / np.abs(z)
    df.loc[df["se"] == 0, "se"] = 1e-10
    return df


# ── Per-locus colocalization ─────────────────────────────────────────────────

def _coloc_one_locus(args):
    """Run coloc for one locus across all ranking cell types."""
    locus, rsid_to_pos, pos_to_rsid, liftover = args
    cfg = GWAS_LOCI[locus]
    chrom = cfg["chr"]
    subtype = cfg["subtype"]
    gene = locus

    # Get genomic position of lead SNP (or proxy)
    lead_snp = cfg["rsid"]
    pos_info = rsid_to_pos.get(lead_snp)
    if pos_info is None:
        for px_id, _ in cfg.get("proxies", []):
            pos_info = rsid_to_pos.get(px_id)
            if pos_info:
                break
    if pos_info is None:
        return [{"locus": locus, "_skip": "lead SNP position unknown"}]

    _, pos_center = pos_info
    pos_center = int(pos_center)

    # Convert lead SNP hg38 position to GWAS build for region query
    gwas_center = pos_center
    if liftover is not None:
        result = liftover.convert_coordinate(f"chr{chrom}", pos_center)
        if result and len(result) > 0:
            # convert_coordinate returns list of (chr, pos, strand, score)
            # but this is hg38->hg38 direction... we need hg38->hg19
            pass
        # Actually we need the reverse: hg38 center -> hg19 for GWAS query
        # Use a reverse liftover
        try:
            from pyliftover import LiftOver
            reverse_lo = LiftOver('hg38', GWAS_BUILD)
            rev_result = reverse_lo.convert_coordinate(f"chr{chrom}", pos_center)
            if rev_result and len(rev_result) > 0:
                gwas_center = int(rev_result[0][1])
        except Exception:
            gwas_center = pos_center  # fallback: positions usually close

    # Load subtype-matched GWAS region (in GWAS build coordinates)
    gwas = load_gwas_region(subtype, chrom, gwas_center)
    if gwas is None or gwas.empty:
        return [{"locus": locus, "_skip": f"no GWAS data for {subtype} at chr{chrom}:{gwas_center}"}]

    # Liftover GWAS positions to hg38 for matching with Bryois
    if liftover is not None:
        hg38_positions = []
        for _, row in gwas.iterrows():
            result = liftover.convert_coordinate(f"chr{row['chr']}", int(row['pos']))
            if result and len(result) > 0:
                hg38_positions.append(int(result[0][1]))
            else:
                hg38_positions.append(None)
        gwas["pos_hg38"] = hg38_positions
        gwas = gwas.dropna(subset=["pos_hg38"])
        gwas["pos_hg38"] = gwas["pos_hg38"].astype(int)
    else:
        gwas["pos_hg38"] = gwas["pos"]

    results = []
    for ct_name, ct_prefix in BRYOIS_CT_PREFIX.items():
        eqtl = load_bryois_region(ct_prefix, chrom, gene, gwas_positions)
        if eqtl is None or eqtl.empty:
            results.append({
                "locus": locus, "gene": gene, "cell_type": ct_name,
                "subtype": subtype, "nsnps": 0,
                "PP.H0": None, "PP.H1": None, "PP.H2": None,
                "PP.H3": None, "PP.H4": None,
                "note": "no eQTL data for gene in this cell type",
            })
            continue

        # Map Bryois rsIDs to positions for merging with GWAS
        eqtl_positions = []
        for rsid in eqtl["snp_id"]:
            p = rsid_to_pos.get(rsid)
            if p:
                eqtl_positions.append((str(p[0]), int(p[1])))
            else:
                eqtl_positions.append((None, None))
        eqtl["chr_pos"] = [f"{c}:{p}" if c else None for c, p in eqtl_positions]
        eqtl["match_chr"] = [c for c, _ in eqtl_positions]
        eqtl["match_pos"] = [p for _, p in eqtl_positions]
        eqtl = eqtl.dropna(subset=["match_pos"])
        eqtl["match_pos"] = eqtl["match_pos"].astype(int)

        # Merge on hg38 position (both sides now in hg38)
        gwas["match_chr"] = gwas["chr"].astype(str)
        gwas["match_pos"] = gwas["pos_hg38"].astype(int)
        merged = eqtl.merge(
            gwas[["match_chr", "match_pos", "beta", "se"]],
            on=["match_chr", "match_pos"],
            suffixes=("_eqtl", "_gwas"),
        )

        if len(merged) < 10:
            results.append({
                "locus": locus, "gene": gene, "cell_type": ct_name,
                "subtype": subtype, "nsnps": len(merged),
                "PP.H0": None, "PP.H1": None, "PP.H2": None,
                "PP.H3": None, "PP.H4": None,
                "note": f"<10 shared SNPs ({len(merged)})",
            })
            continue

        # Drop rows with zero/nan SE
        merged = merged[
            (merged["se_gwas"] > 0) & (merged["se_eqtl"] > 0) &
            merged["se_gwas"].notna() & merged["se_eqtl"].notna()
        ]

        if len(merged) < 10:
            results.append({
                "locus": locus, "gene": gene, "cell_type": ct_name,
                "subtype": subtype, "nsnps": len(merged),
                "PP.H0": None, "PP.H1": None, "PP.H2": None,
                "PP.H3": None, "PP.H4": None,
                "note": f"<10 valid SNPs after SE filter ({len(merged)})",
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
        pp["subtype"] = subtype
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
    missing = []
    for subtype, fname in GWAS_SUBTYPE_FILES.items():
        if not (GWAS_DIR / fname).exists():
            missing.append(f"  {subtype}: {GWAS_DIR / fname}")

    if missing:
        log.info("=" * 72)
        log.info("GWAS summary statistics not found — skipping colocalization")
        log.info("=" * 72)
        log.info("")
        log.info("Missing files:")
        for m in missing:
            log.info(m)
        log.info("")
        log.info("Create symlinks in data/gwas/ to your meta-analysis files:")
        log.info("  cd data/gwas")
        for subtype, fname in GWAS_SUBTYPE_FILES.items():
            log.info(f"  ln -s /path/to/{fname} {fname}")
        log.info("")
        log.info("The pipeline continues without colocalization.")

        stub = pd.DataFrame(columns=["locus", "gene", "cell_type",
                                      "PP.H4", "note"])
        stub.to_csv(OUTPUT_DIR / "coloc_results.csv", index=False)
        return

    log.info("=" * 72)
    log.info("Colocalization analysis (coloc.abf)")
    log.info("=" * 72)
    for subtype, fname in GWAS_SUBTYPE_FILES.items():
        fpath = GWAS_DIR / fname
        size_mb = fpath.stat().st_size / 1e6
        log.info(f"  {subtype:10s}: {fname} ({size_mb:.0f} MB)")
    log.info(f"  Using {ncpus} workers")

    # Load SNP positions for coordinate-based matching
    log.info("  Loading Bryois SNP positions (hg38)...")
    rsid_to_pos, pos_to_rsid = _get_snp_positions()
    log.info(f"  {len(rsid_to_pos):,} rsID->position mappings loaded")

    # Set up liftover if GWAS is on a different build
    liftover = None
    if GWAS_BUILD != "hg38":
        try:
            from pyliftover import LiftOver
            liftover = LiftOver(GWAS_BUILD, 'hg38')
            log.info(f"  Liftover: {GWAS_BUILD} -> hg38 (via pyliftover)")
        except ImportError:
            log.error("  pyliftover not installed! pip install pyliftover")
            log.error("  Cannot match GWAS (hg19) to Bryois (hg38) without liftover.")
            return
    else:
        log.info("  GWAS is hg38 — no liftover needed")

    # Run per-locus in parallel
    tasks = [(locus, rsid_to_pos, pos_to_rsid, liftover) for locus in GWAS_LOCI]
    all_results = []
    with ProcessPoolExecutor(max_workers=ncpus) as pool:
        futures = {pool.submit(_coloc_one_locus, t): t[0] for t in tasks}
        for fut in as_completed(futures):
            locus = futures[fut]
            try:
                rows = fut.result()
            except Exception as e:
                log.error(f"  [{locus}] error: {e}")
                continue
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
        summary = valid.loc[idx, ["locus", "gene", "cell_type", "subtype",
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
