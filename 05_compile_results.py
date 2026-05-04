#!/usr/bin/env python3
"""
05_compile_results.py
---------------------
Merge all upstream outputs into a comprehensive per-locus table and write
an interpretive markdown report.

Output
------
  output/comprehensive_results.csv   – one row per locus, all key columns
  output/interpretation.md           – markdown report with project context
"""

import numpy as np
import pandas as pd
from datetime import datetime

from config import (
    GWAS_LOCI, OUTPUT_DIR, RANKING_CELL_TYPES,
    BRYOIS_CT_PREFIX, NOMINAL_P, BH_FDR, COLOC_PP4,
    setup_logging,
)

log = setup_logging("05-compile", OUTPUT_DIR / "05_compile.log")


def build_comprehensive_table():
    log.info("=" * 72)
    log.info("Building comprehensive results table")
    log.info("=" * 72)

    # Load upstream results (all optional — degrade gracefully)
    expr_df = _load("celltype_expression.csv")
    summary_df = _load("locus_summary.csv")
    eqtl_df = _load("eqtl_corrected.csv")
    coloc_df = _load("coloc_summary.csv")

    rows = []
    for locus, cfg in GWAS_LOCI.items():
        row = {
            "locus":       locus,
            "rsid":        cfg["rsid"],
            "chr":         cfg["chr"],
            "locus_band":  cfg["locus"],
            "subtype":     cfg["subtype"],
            "pathway":     cfg["pathway"],
        }

        # Expression in OPC and Oligodendrocyte
        if expr_df is not None:
            for ct in ["OPC", "Oligodendrocyte"]:
                r = expr_df[(expr_df["gene"] == locus) &
                            (expr_df["cell_type"] == ct)]
                col = f"pct_expr_{ct.lower()}"
                row[col] = round(r["pct_expressing"].iloc[0], 1) if len(r) else None

        # Best eQTL from summary
        if summary_df is not None:
            s = summary_df[summary_df["locus"] == locus]
            if len(s):
                r = s.iloc[0]
                row["best_celltype"]     = r.get("best_cell_type")
                row["best_pvalue"]       = r.get("best_pvalue")
                row["best_pvalue_bh"]    = r.get("best_pvalue_bh")
                row["best_beta"]         = r.get("best_beta")
                row["snp_used"]          = r.get("snp_used")
                row["survives_bh"]       = r.get("best_survives_bh")
                row["n_sig_nominal"]     = r.get("n_sig_nominal")
                row["n_sig_bh"]          = r.get("n_sig_bh")
                if pd.notna(r.get("concordant_direction")):
                    row["concordant_direction"] = r["concordant_direction"]

        # Per-cell-type p-values (BH-corrected)
        if eqtl_df is not None:
            for ct in RANKING_CELL_TYPES:
                r = eqtl_df[(eqtl_df["locus"] == locus) &
                            (eqtl_df["cell_type"] == ct)]
                if len(r):
                    rb = r.loc[r["pvalue"].idxmin()]
                    row[f"p_{ct}"] = rb["pvalue"]
                    if "pvalue_bh" in rb.index:
                        row[f"p_bh_{ct}"] = rb["pvalue_bh"]

        # Colocalization
        if coloc_df is not None:
            c = coloc_df[coloc_df["locus"] == locus]
            if len(c):
                row["coloc_best_ct"]  = c.iloc[0].get("cell_type")
                row["coloc_PP_H4"]    = c.iloc[0].get("PP.H4")
                row["colocalized"]    = c.iloc[0].get("colocalized")

        rows.append(row)

    comp = pd.DataFrame(rows)
    out = OUTPUT_DIR / "comprehensive_results.csv"
    comp.to_csv(out, index=False)
    log.info(f"  Saved {out} ({len(comp)} loci)")
    return comp


def write_interpretation(comp_df):
    log.info("=" * 72)
    log.info("Writing interpretation report")
    log.info("=" * 72)

    summary_df = _load("locus_summary.csv")
    perm_df = _load("permutation_results.csv")
    coloc_summary = _load("coloc_summary.csv")

    lines = [
        "# Single-Cell eQTL Analysis: Locus-Level Results",
        "",
        f"_Generated: {datetime.now().isoformat(timespec='seconds')}_",
        "",
        "## Pipeline improvements over initial script",
        "",
        "1. **Within-locus BH correction** across cell types (FDR < 0.05)",
        "2. **Single best proxy** per locus (highest r²) — no cherry-picking",
        "3. **Pseudobulk excluded** from cell-type ranking (kept as reference)",
        "4. **Risk allele direction** assessed where GWAS betas available",
        "5. **Colocalization** (coloc.abf) tests shared causal variants",
        "6. **Permutation test** for subtype × cell-type enrichment",
        "",
        "## How to read these results",
        "",
        "This locus-level analysis is the **companion** to the cell-stratified",
        "MR (csMR). The csMR tests ICVF → glioma causality at the *trait level*.",
        "This script tests each individual GWAS lead SNP to identify which",
        "cell type its regulatory effect is strongest in.",
        "",
        "**The headline csMR finding (replicated in FinnGen):**",
        "- Oligodendrocyte-regulatory ICVF variants → IDH-wt glioma",
        "- OR=1.82, p=0.004 (FinnGen: OR=2.01, p=0.001)",
        "",
        "Results here are **supporting examples**, not independent evidence",
        "with comparable rigor.",
        "",
        "## Locus-level summary",
        "",
    ]

    if summary_df is not None:
        lines += [
            "| Locus | Subtype | Best Cell Type | p (raw) | p (BH) | β | "
            "Survives BH | SNP |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for _, r in summary_df.sort_values(["subtype", "best_pvalue"]).iterrows():
            if pd.isna(r.get("best_pvalue")):
                lines.append(
                    f"| {r['locus']} | {r['subtype']} | (not testable) "
                    f"| - | - | - | - | - |")
            else:
                surv = "✓" if r.get("best_survives_bh") else "✗"
                p_bh = (f"{r['best_pvalue_bh']:.2e}"
                        if pd.notna(r.get("best_pvalue_bh")) else "N/A")
                lines.append(
                    f"| {r['locus']} | {r['subtype']} | "
                    f"{r['best_cell_type']} | {r['best_pvalue']:.2e} | "
                    f"{p_bh} | "
                    f"{r['best_beta']:+.3f} | {surv} | {r['snp_used']} |")

    # Colocalization section
    lines += ["", "## Colocalization results", ""]
    if coloc_summary is not None and not coloc_summary.empty:
        lines += [
            "| Locus | Best Cell Type | PP.H4 | Colocalized |",
            "|---|---|---|---|",
        ]
        for _, r in coloc_summary.iterrows():
            flag = "✓" if r.get("colocalized") else "✗"
            lines.append(
                f"| {r['locus']} | {r.get('cell_type', '-')} | "
                f"{r.get('PP.H4', 0):.3f} | {flag} |")
    else:
        lines.append("_Colocalization not run (GWAS summary stats not provided)._")

    # Permutation results
    lines += ["", "## Enrichment permutation tests", ""]
    if perm_df is not None and not perm_df.empty:
        lines += ["| Test | Focus prop | Other prop | Δ | Perm p |",
                   "|---|---|---|---|---|"]
        for _, r in perm_df.iterrows():
            lines.append(
                f"| {r['test']} | {r['obs_prop_focus']:.3f} | "
                f"{r['obs_prop_other']:.3f} | {r['obs_diff']:+.3f} | "
                f"{r['perm_p']:.4f} |")
        lines.append("")
        lines.append(f"_Based on {int(perm_df['n_perm'].iloc[0]):,} permutations. "
                     f"n=12 loci — interpret as descriptive._")
    else:
        lines.append("_Permutation results not available._")

    # Caveats
    lines += [
        "", "## Caveats", "",
        "1. **N=12 loci** is too small for formal group-level tests.",
        "2. **LD proxies** for CCDC26 (max r²=0.64) and CDKN2A "
        "(r² unknown) attenuate effect sizes.",
        "3. **Bryois N=196** limits power for less-abundant cell types.",
        "4. **No replication** at the single-cell eQTL level.",
        "5. **Within-locus BH correction** is conservative for n=6 tests "
        "— some real effects may be FDR-corrected away.",
        "",
        "## Output files",
        "",
        "| File | Description |",
        "|---|---|",
        "| `celltype_expression.csv` | % expressing + mean, 12 genes × 8 types |",
        "| `eqtl_results_all.csv` | Long-format: every (locus, snp, gene, ct) |",
        "| `eqtl_corrected.csv` | BH-corrected within-locus results |",
        "| `locus_summary.csv` | Best cell type per locus (corrected) |",
        "| `permutation_results.csv` | Enrichment permutation tests |",
        "| `coloc_results.csv` | Full coloc output (if GWAS data provided) |",
        "| `coloc_summary.csv` | Best PP.H4 per locus |",
        "| `comprehensive_results.csv` | One-row-per-locus master table |",
        "| `fig1–fig4` | Publication-ready figures |",
        "",
        "## Reference",
        "",
        "Bryois J et al. Cell-type-specific cis-eQTLs in eight human brain",
        "cell types identify novel risk genes for psychiatric and neurological",
        "disorders. *Nat Neurosci* 25, 1104–1112 (2022).",
        "DOI: 10.1038/s41593-022-01128-z",
    ]

    out = OUTPUT_DIR / "interpretation.md"
    out.write_text("\n".join(lines))
    log.info(f"  Saved {out}")


def _load(name):
    path = OUTPUT_DIR / name
    if path.exists():
        return pd.read_csv(path)
    return None


def main():
    comp = build_comprehensive_table()
    write_interpretation(comp)
    log.info("Done.")


if __name__ == "__main__":
    main()
