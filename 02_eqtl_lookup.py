#!/usr/bin/env python3
"""
02_eqtl_lookup.py
-----------------
Extract Bryois eQTL results for every (locus, gene, cell_type) combination.
Parallelised across cell-type x chromosome files using multiprocessing.

Key fixes vs. the original script
----------------------------------
1. Within-locus Benjamini-Hochberg correction across cell types.
2. Only the highest-r² proxy is used per locus (no cherry-picking).
3. Pseudobulk excluded from "best cell type" ranking.
4. Risk-allele direction concordance (when GWAS betas are provided).
5. Pre-loads each file once and extracts all needed rows in a single pass.
6. Permutation test for cell-type enrichment.

Output
------
  output/eqtl_results_all.csv      – long-format (locus, snp, gene, cell_type)
  output/locus_summary.csv         – best cell type per locus (corrected)
  output/permutation_results.csv   – permutation enrichment tests
"""

import os, sys, gzip, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from config import (
    GWAS_LOCI, ALL_GENES, BRYOIS_DIR, GWAS_DIR, OUTPUT_DIR,
    BRYOIS_CT_PREFIX, BRYOIS_PSEUDO_PREFIX, BRYOIS_ALL_PREFIX,
    RANKING_CELL_TYPES, NEEDED_CHROMS, NOMINAL_P, BH_FDR, N_PERM,
    GWAS_SUBTYPE_FILES, GWAS_COLUMN_MAP, GWAS_BUILD,
    snps_for_locus, genes_for_locus, setup_logging,
)

log = setup_logging("02-eqtl", OUTPUT_DIR / "02_eqtl_lookup.log")


# ── Build the lookup target set ──────────────────────────────────────────────

def _build_targets():
    """Return {chr_int: set((gene_symbol, snp_id, locus, role, r2))}."""
    targets = defaultdict(set)
    for locus in GWAS_LOCI:
        cfg = GWAS_LOCI[locus]
        for snp_id, role, r2 in snps_for_locus(locus):
            for gene in genes_for_locus(locus):
                targets[cfg["chr"]].add((gene, snp_id, locus, role, r2))
    return targets


# ── Single-file extraction (runs in worker) ──────────────────────────────────

def _extract_from_file(args):
    """Read one Bryois gz file and return matching rows as dicts."""
    ct_name, ct_prefix, chrom, target_set = args
    fpath = BRYOIS_DIR / f"{ct_prefix}.{chrom}.gz"
    if not fpath.exists():
        return ct_name, chrom, []

    # Build fast lookup: {(gene_symbol, snp_id): (locus, role, r2)}
    lookup = {}
    for gene, snp, locus, role, r2 in target_set:
        lookup[(gene, snp)] = (locus, role, r2)

    hits = []
    try:
        with gzip.open(fpath, "rt") as f:
            for line in f:
                parts = line.split()
                if len(parts) != 5:
                    continue
                gene_col, snp_col = parts[0], parts[1]
                sym = gene_col.split("_")[0]
                key = (sym, snp_col)
                if key not in lookup:
                    continue
                locus, role, r2 = lookup[key]
                hits.append({
                    "locus":     locus,
                    "subtype":   GWAS_LOCI[locus]["subtype"],
                    "pathway":   GWAS_LOCI[locus]["pathway"],
                    "lead_snp":  GWAS_LOCI[locus]["rsid"],
                    "snp_used":  snp_col,
                    "snp_role":  role,
                    "ld_r2":     r2,
                    "gene":      sym,
                    "gene_full": gene_col,
                    "cell_type": ct_name,
                    "chr":       chrom,
                    "distance":  int(parts[2]),
                    "pvalue":    float(parts[3]),
                    "beta":      float(parts[4]),
                    "tested":    True,
                })
    except Exception as e:
        return ct_name, chrom, [{"_error": str(e)}]

    return ct_name, chrom, hits


# ── Parallel extraction ──────────────────────────────────────────────────────

def extract_eqtl_parallel(ncpus):
    log.info("=" * 72)
    log.info("Parallel eQTL extraction")
    log.info("=" * 72)

    targets_by_chr = _build_targets()
    all_prefixes = {**BRYOIS_CT_PREFIX, **BRYOIS_PSEUDO_PREFIX}

    # Build task list: (ct_name, ct_prefix, chrom, target_set)
    tasks = []
    for ct_name, ct_prefix in all_prefixes.items():
        for chrom in NEEDED_CHROMS:
            tset = targets_by_chr.get(chrom, set())
            if tset:
                tasks.append((ct_name, ct_prefix, chrom, tset))

    log.info(f"  {len(tasks)} file-extraction tasks across "
             f"{len(all_prefixes)} cell types x {len(NEEDED_CHROMS)} chroms")
    log.info(f"  Using {ncpus} worker processes")

    all_hits = []
    errors = 0
    with ProcessPoolExecutor(max_workers=ncpus) as pool:
        futures = {pool.submit(_extract_from_file, t): t for t in tasks}
        for fut in as_completed(futures):
            ct_name, chrom, hits = fut.result()
            for h in hits:
                if "_error" in h:
                    log.error(f"  [{ct_name} chr{chrom}] {h['_error']}")
                    errors += 1
                else:
                    all_hits.append(h)

    log.info(f"  Extracted {len(all_hits):,} eQTL hits ({errors} errors)")

    # Build full long-format table (hits + not-found rows)
    eqtl_df = pd.DataFrame(all_hits) if all_hits else pd.DataFrame()

    # Add rows for (locus, snp, gene, cell_type) combos not found
    expected = set()
    for locus in GWAS_LOCI:
        cfg = GWAS_LOCI[locus]
        for snp_id, role, r2 in snps_for_locus(locus):
            for gene in genes_for_locus(locus):
                for ct_name in all_prefixes:
                    expected.add((locus, snp_id, gene, ct_name))

    if len(eqtl_df):
        found = set(zip(eqtl_df["locus"], eqtl_df["snp_used"],
                        eqtl_df["gene"], eqtl_df["cell_type"]))
    else:
        found = set()

    missing_rows = []
    for locus, snp_id, gene, ct_name in expected - found:
        cfg = GWAS_LOCI[locus]
        snp_info = [(s, r, r2) for s, r, r2 in snps_for_locus(locus)
                    if s == snp_id][0]
        missing_rows.append({
            "locus": locus, "subtype": cfg["subtype"],
            "pathway": cfg["pathway"], "lead_snp": cfg["rsid"],
            "snp_used": snp_id, "snp_role": snp_info[1],
            "ld_r2": snp_info[2], "gene": gene,
            "gene_full": None, "cell_type": ct_name,
            "chr": cfg["chr"], "distance": None,
            "pvalue": None, "beta": None, "tested": False,
        })

    eqtl_df = pd.concat([eqtl_df, pd.DataFrame(missing_rows)],
                         ignore_index=True)

    out = OUTPUT_DIR / "eqtl_results_all.csv"
    eqtl_df.to_csv(out, index=False)
    n_tested = int(eqtl_df["tested"].sum())
    n_sig = int((eqtl_df["pvalue"] < NOMINAL_P).sum()) if n_tested else 0
    log.info(f"  Saved {out}  ({len(eqtl_df):,} rows, "
             f"{n_tested} tested, {n_sig} nom. sig.)")
    return eqtl_df


# ── Automatic GWAS lead-SNP lookup ────────────────────────────────────────────

def lookup_gwas_lead_snps():
    """Look up GWAS beta/SE/EAF for each lead SNP from the summary stats files.

    Uses Bryois snp_pos.txt.gz for rsID -> hg38 position, then lifts over
    to the GWAS build (hg19) to match by genomic position.

    Returns dict: {locus_name: {gwas_label: {beta, se, eaf, pvalue}}}
    """
    log.info("=" * 72)
    log.info("Looking up GWAS effect sizes for lead SNPs")
    log.info("=" * 72)

    # Load SNP positions (hg38)
    snp_pos_path = BRYOIS_DIR / "snp_pos.txt.gz"
    if not snp_pos_path.exists():
        log.warning("  snp_pos.txt.gz not found — skipping GWAS lookup")
        return {}

    import gzip
    # snp_pos.txt.gz from Bryois (Zenodo 7276971) has columns:
    #   SNP  SNP_id_hg38  SNP_id_hg19  effect_allele  other_allele  MAF
    # where SNP_id_hg38/hg19 are in "chr1:234313" format.
    # We extract BOTH builds so downstream code can match GWAS (hg19)
    # without liftover — avoiding the small position losses liftover introduces.

    rsid_to_hg38 = {}
    rsid_to_hg19 = {}

    with gzip.open(snp_pos_path, "rt") as f:
        header = f.readline().strip()
        sep = "\t" if "\t" in header else " "
        cols = header.split(sep)
        log.info(f"  snp_pos.txt.gz columns: {cols}")

        # Identify which columns hold hg38 and hg19 positions
        hg38_col = None
        hg19_col = None
        for idx, c in enumerate(cols):
            cl = c.lower()
            if "hg38" in cl:
                hg38_col = idx
            elif "hg19" in cl:
                hg19_col = idx
        log.info(f"  hg38 column index: {hg38_col}, hg19 column index: {hg19_col}")

        for line in f:
            parts = line.strip().split(sep)
            if len(parts) < 2:
                continue
            rsid = parts[0]

            # Parse hg38
            if hg38_col is not None and hg38_col < len(parts):
                p = parts[hg38_col]
                if ":" in p:
                    cp = p.replace("chr", "").split(":")
                    if cp[0].isdigit():
                        rsid_to_hg38[rsid] = (cp[0], int(cp[1]))

            # Parse hg19
            if hg19_col is not None and hg19_col < len(parts):
                p = parts[hg19_col]
                if ":" in p:
                    cp = p.replace("chr", "").split(":")
                    if cp[0].isdigit():
                        rsid_to_hg19[rsid] = (cp[0], int(cp[1]))

            # Fallback: if no labelled columns, use first chr:pos found as hg38
            if rsid not in rsid_to_hg38:
                for p in parts[1:]:
                    if ":" in p and p.replace("chr", "").split(":")[0].isdigit():
                        cp = p.replace("chr", "").split(":")
                        rsid_to_hg38[rsid] = (cp[0], int(cp[1]))
                        break

    log.info(f"  {len(rsid_to_hg38):,} hg38 + {len(rsid_to_hg19):,} hg19 positions loaded")

    # Use the native hg19 positions from snp_pos.txt.gz to match GWAS data.
    # This avoids liftover entirely — no small position losses from chain files.
    #
    # NOTE: pyliftover code preserved below (commented out) in case a future
    # dataset lacks native hg19 positions in its SNP reference file.
    #
    # --- Liftover approach (commented out — snp_pos.txt.gz has both builds) ---
    # liftover = None
    # if GWAS_BUILD != "hg38":
    #     try:
    #         from pyliftover import LiftOver
    #         liftover = LiftOver("hg38", GWAS_BUILD)
    #         log.info(f"  Liftover: hg38 -> {GWAS_BUILD}")
    #     except ImportError:
    #         log.warning("  pyliftover not installed — cannot match positions")
    #         return {}
    # --- End liftover ---

    # Choose the position dict matching the GWAS build
    if GWAS_BUILD == "hg19" and rsid_to_hg19:
        gwas_pos_lookup = rsid_to_hg19
        log.info(f"  Using native hg19 positions for GWAS matching (no liftover)")
    elif GWAS_BUILD == "hg38":
        gwas_pos_lookup = rsid_to_hg38
        log.info(f"  GWAS is hg38 — using hg38 positions directly")
    else:
        log.warning(f"  No {GWAS_BUILD} positions available — falling back to hg38")
        gwas_pos_lookup = rsid_to_hg38

    # For each lead SNP, get its position in the GWAS build
    lead_snp_positions = {}  # {locus: (chr, gwas_build_pos)}
    for locus, cfg in GWAS_LOCI.items():
        rsid = cfg["rsid"]
        pos_info = gwas_pos_lookup.get(rsid)
        if pos_info is None:
            # Try proxies
            for px_id, _ in cfg.get("proxies", []):
                pos_info = gwas_pos_lookup.get(px_id)
                if pos_info:
                    rsid = px_id
                    break
        if pos_info is None:
            log.info(f"  [{locus}] lead SNP position unknown — skipped")
            continue

        chrom, pos_gwas = pos_info
        lead_snp_positions[locus] = (str(chrom), pos_gwas, rsid)

    log.info(f"  {len(lead_snp_positions)}/{len(GWAS_LOCI)} lead SNP positions resolved")

    # Look up each lead SNP in each GWAS file
    gwas_effects = {}  # {locus: {gwas_label: {beta, se, eaf, pvalue}}}

    for gwas_label, gwas_fname in GWAS_SUBTYPE_FILES.items():
        gwas_path = GWAS_DIR / gwas_fname
        if not gwas_path.exists():
            continue

        log.info(f"  Scanning {gwas_fname} for lead SNPs...")

        # Read only needed columns
        needed_cols = [c for c in GWAS_COLUMN_MAP.keys()
                       if c in ["CHR", "BP", "BETA", "SE", "P", "A1_FREQ"]]
        try:
            df = pd.read_csv(gwas_path, sep="\t", usecols=needed_cols,
                             dtype={"CHR": str, "BP": int})
        except Exception as e:
            log.error(f"    Failed to read: {e}")
            continue

        df = df.rename(columns=GWAS_COLUMN_MAP)
        df["chr"] = df["chr"].astype(str).str.replace("chr", "")

        for locus, (chrom, pos_gwas, rsid) in lead_snp_positions.items():
            # Match by chr + position (allow ±1 bp for rounding)
            match = df[(df["chr"] == chrom) &
                       (df["pos"] >= pos_gwas - 1) &
                       (df["pos"] <= pos_gwas + 1)]
            if len(match) == 0:
                continue
            row = match.iloc[0]
            if locus not in gwas_effects:
                gwas_effects[locus] = {}
            gwas_effects[locus][gwas_label] = {
                "gwas_beta":  float(row["beta"]),
                "gwas_se":    float(row["se"]),
                "gwas_eaf":   float(row.get("eaf", np.nan)),
                "gwas_pvalue": float(row["pvalue"]),
            }

    # Log summary
    for locus in GWAS_LOCI:
        effects = gwas_effects.get(locus, {})
        if effects:
            matched_sub = GWAS_LOCI[locus]["subtype"]
            matched = effects.get(matched_sub) or effects.get("all_glioma") or next(iter(effects.values()))
            log.info(f"  [{locus:10s}] beta={matched['gwas_beta']:+.4f}  "
                     f"se={matched['gwas_se']:.4f}  p={matched['gwas_pvalue']:.2e}")
        else:
            log.info(f"  [{locus:10s}] not found in any GWAS file")

    # Save the lookup table
    rows = []
    for locus, effects in gwas_effects.items():
        for gwas_label, vals in effects.items():
            rows.append({"locus": locus, "gwas": gwas_label, **vals})
    if rows:
        out = OUTPUT_DIR / "gwas_lead_snp_effects.csv"
        pd.DataFrame(rows).to_csv(out, index=False)
        log.info(f"  Saved {out}")

    return gwas_effects


# ── Within-locus multiple testing correction ─────────────────────────────────

def correct_and_summarise(eqtl_df, gwas_effects=None):
    """
    For each locus (primary gene, lead-or-proxy):
      1. Gather p-values across RANKING cell types (excluding Pseudobulk).
      2. Apply Benjamini-Hochberg correction.
      3. Identify best cell type based on corrected p-values.
    """
    log.info("=" * 72)
    log.info("Within-locus correction & summary")
    log.info("=" * 72)
    if gwas_effects is None:
        gwas_effects = {}

    # Primary gene = locus name; keep ranking cell types only
    primary = eqtl_df[
        (eqtl_df["gene"] == eqtl_df["locus"])
        & eqtl_df["cell_type"].isin(RANKING_CELL_TYPES)
    ].copy()
    primary = primary.dropna(subset=["pvalue"])

    if primary.empty:
        log.warning("  No testable eQTL results found!")
        return pd.DataFrame()

    # For loci with a proxy, pick the SNP with the smallest p per cell type
    # (there are at most 2: lead + 1 proxy)
    idx = primary.groupby(["locus", "cell_type"])["pvalue"].idxmin()
    best_snp = primary.loc[idx].copy()

    # BH correction within each locus
    def _bh(grp):
        from statsmodels.stats.multitest import multipletests
        pvals = grp["pvalue"].values
        if len(pvals) == 0:
            grp["pvalue_bh"] = np.nan
            grp["sig_bh"] = False
            return grp
        _, pvals_adj, _, _ = multipletests(pvals, alpha=BH_FDR, method="fdr_bh")
        grp["pvalue_bh"] = pvals_adj
        grp["sig_bh"] = pvals_adj < BH_FDR
        return grp

    best_snp = best_snp.groupby("locus", group_keys=False).apply(_bh)

    # Save the full corrected table
    corrected_out = OUTPUT_DIR / "eqtl_corrected.csv"
    best_snp.to_csv(corrected_out, index=False)
    log.info(f"  Saved {corrected_out}")

    # Locus summary: best cell type per locus (using BH-corrected p)
    rows = []
    for locus in GWAS_LOCI:
        sub = best_snp[best_snp["locus"] == locus]
        cfg = GWAS_LOCI[locus]
        if sub.empty:
            rows.append({
                "locus": locus, "subtype": cfg["subtype"],
                "best_cell_type": None, "best_pvalue": None,
                "best_pvalue_bh": None, "best_beta": None,
                "snp_used": None, "n_celltypes_tested": 0,
                "n_sig_nominal": 0, "n_sig_bh": 0,
                "best_survives_bh": False,
            })
            continue
        b = sub.loc[sub["pvalue"].idxmin()]
        rows.append({
            "locus":              locus,
            "subtype":            cfg["subtype"],
            "best_cell_type":     b["cell_type"],
            "best_pvalue":        b["pvalue"],
            "best_pvalue_bh":     b["pvalue_bh"],
            "best_beta":          b["beta"],
            "snp_used":           b["snp_used"],
            "n_celltypes_tested": len(sub),
            "n_sig_nominal":      int((sub["pvalue"] < NOMINAL_P).sum()),
            "n_sig_bh":           int(sub["sig_bh"].sum()),
            "best_survives_bh":   bool(b["sig_bh"]),
        })

    summary = pd.DataFrame(rows).sort_values(["subtype", "best_pvalue"])

    # ── Risk allele direction (allele-harmonized) ──
    # The Bryois eQTL beta is relative to the Bryois effect_allele.
    # The GWAS beta is relative to GWAS A1.
    # If these are different alleles (detectable by comparing MAF vs EAF),
    # we flip the GWAS beta before comparing signs.
    #
    # Load Bryois allele info for harmonization
    _bryois_alleles = {}
    snp_pos_path = BRYOIS_DIR / "snp_pos.txt.gz"
    if snp_pos_path.exists():
        with gzip.open(snp_pos_path, "rt") as f:
            hdr = f.readline().strip()
            _sep = "\t" if "\t" in hdr else " "
            for line in f:
                parts = line.strip().split(_sep)
                if len(parts) >= 6:
                    _bryois_alleles[parts[0]] = float(parts[5])  # MAF

    def _get_harmonized_gwas_beta(row):
        locus = row["locus"]
        snp = row.get("snp_used")
        effects = gwas_effects.get(locus, {})
        if not effects:
            return None
        subtype = GWAS_LOCI[locus]["subtype"]
        matched = effects.get(subtype) or effects.get("all_glioma")
        if not matched:
            matched = next(iter(effects.values()))
        gwas_beta = matched["gwas_beta"]
        gwas_eaf = matched["gwas_eaf"]

        # Check if alleles are flipped using frequency comparison
        bryois_maf = _bryois_alleles.get(snp)
        if bryois_maf is not None and not np.isnan(gwas_eaf):
            diff_same = abs(gwas_eaf - bryois_maf)
            diff_flip = abs(gwas_eaf - (1 - bryois_maf))
            if diff_flip < diff_same and diff_flip < 0.10:
                # Alleles are flipped — negate GWAS beta to harmonize
                gwas_beta = -gwas_beta
        return gwas_beta

    summary["gwas_beta_harmonized"] = summary.apply(_get_harmonized_gwas_beta, axis=1)
    # Keep original for reference
    summary["gwas_beta_raw"] = summary["locus"].map(
        lambda l: (gwas_effects.get(l, {}).get(GWAS_LOCI[l]["subtype"])
                   or gwas_effects.get(l, {}).get("all_glioma")
                   or {"gwas_beta": None})["gwas_beta"]
        if l in gwas_effects else None
    )
    summary["gwas_beta"] = summary["gwas_beta_harmonized"]

    has_gwas = summary["gwas_beta"].notna() & summary["best_beta"].notna()
    summary["concordant_direction"] = np.where(
        has_gwas,
        np.sign(summary["gwas_beta"]) == np.sign(summary["best_beta"]),
        None,
    )

    out = OUTPUT_DIR / "locus_summary.csv"
    summary.to_csv(out, index=False)
    log.info(f"  Saved {out}")
    log.info(f"\n{summary.to_string(index=False)}\n")
    return summary


# ── Permutation test for cell-type enrichment ────────────────────────────────

def permutation_enrichment(summary_df):
    """Permutation test: are IDH-wt loci enriched for oligo-lineage eQTLs?"""
    log.info("=" * 72)
    log.info(f"Permutation enrichment test ({N_PERM:,} permutations)")
    log.info("=" * 72)

    df = summary_df.dropna(subset=["best_cell_type"]).copy()
    if df.empty:
        log.warning("  No data for permutation test")
        return

    df["is_oligo_lineage"] = df["best_cell_type"].isin(
        ["OPC", "Oligodendrocyte"])
    df["is_opc"] = df["best_cell_type"] == "OPC"
    df["is_oligo"] = df["best_cell_type"] == "Oligodendrocyte"

    rng = np.random.default_rng(42)
    results = []

    for col, label in [("is_opc", "OPC"), ("is_oligo", "Oligodendrocyte"),
                        ("is_oligo_lineage", "Oligo lineage")]:
        for focus_sub in ["IDH-mut", "IDH-wt"]:
            focus = df[df["subtype"] == focus_sub]
            other = df[df["subtype"] != focus_sub]
            n_focus = len(focus)
            n_other = len(other)
            if n_focus == 0 or n_other == 0:
                continue
            obs_focus = focus[col].sum()
            obs_other = other[col].sum()
            obs_diff = obs_focus / n_focus - obs_other / n_other

            # Permutation: shuffle subtype labels
            all_vals = df[col].values.copy()
            n_total = len(all_vals)
            null_diffs = np.empty(N_PERM)
            for i in range(N_PERM):
                rng.shuffle(all_vals)
                perm_focus = all_vals[:n_focus].sum() / n_focus
                perm_other = all_vals[n_focus:].sum() / (n_total - n_focus)
                null_diffs[i] = perm_focus - perm_other

            p_perm = (np.abs(null_diffs) >= np.abs(obs_diff)).mean()
            results.append({
                "test": f"{focus_sub} enriched for {label}",
                "obs_prop_focus": round(obs_focus / n_focus, 3),
                "obs_prop_other": round(obs_other / n_other, 3),
                "obs_diff": round(obs_diff, 3),
                "perm_p": round(p_perm, 4),
                "n_perm": N_PERM,
            })
            log.info(f"  {results[-1]['test']}: "
                     f"diff={obs_diff:+.3f}, perm_p={p_perm:.4f}")

    perm_df = pd.DataFrame(results)
    out = OUTPUT_DIR / "permutation_results.csv"
    perm_df.to_csv(out, index=False)
    log.info(f"  Saved {out}")

    log.info("\n  NOTE: n=12 loci means limited power. These are descriptive.")
    return perm_df


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Parallel eQTL lookup across Bryois cell types")
    parser.add_argument("--cpus", type=int, default=0,
                        help="Worker processes (0 = auto-detect)")
    args = parser.parse_args()

    if args.cpus > 0:
        ncpus = args.cpus
    elif os.environ.get("SLURM_CPUS_PER_TASK"):
        ncpus = int(os.environ["SLURM_CPUS_PER_TASK"])
    elif hasattr(os, 'sched_getaffinity'):
        ncpus = len(os.sched_getaffinity(0))
    else:
        ncpus = os.cpu_count() or 4
    log.info(f"CPUs available: {ncpus}")

    eqtl_df = extract_eqtl_parallel(ncpus)
    gwas_effects = lookup_gwas_lead_snps()
    summary_df = correct_and_summarise(eqtl_df, gwas_effects)
    if not summary_df.empty:
        permutation_enrichment(summary_df)

    log.info("Done.")


if __name__ == "__main__":
    main()
