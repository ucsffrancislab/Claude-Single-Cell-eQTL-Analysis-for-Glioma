#!/usr/bin/env python3
"""
04_visualization.py
-------------------
Generate all figures from the pipeline outputs.

Figures
-------
  fig1_expression_dotplot.png   – 12 genes x 8 cell types (split by subtype)
  fig2_eqtl_heatmap.png        – eQTL signed -log10(p) heatmap (BH-corrected)
  fig3_locus_summary.png        – best cell type per locus, horizontal bars
  fig4_coloc_heatmap.png        – PP.H4 heatmap (only if coloc ran)

All plots use fig.savefig (not plt.savefig) for correct lineage tracking.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib import cm

from config import (
    GWAS_LOCI, OUTPUT_DIR,
    EXPRESSION_CELL_TYPES, RANKING_CELL_TYPES,
    BRYOIS_CT_PREFIX, COLOC_PP4,
    setup_logging,
)

log = setup_logging("04-viz", OUTPUT_DIR / "04_visualization.log")

COLOR_MUT = "#d62728"
COLOR_WT  = "#1f77b4"

CT_COLORS = {
    "OPC":              "#ff7f00",
    "Oligodendrocyte":  "#e31a1c",
    "Astrocyte":        "#33a02c",
    "ExcitatoryNeuron": "#1f78b4",
    "InhibitoryNeuron": "#a6cee3",
    "Microglia":        "#6a3d9a",
}


def _subtype_color(s):
    return COLOR_MUT if s == "IDH-mut" else COLOR_WT


# ── Fig 1: Expression dot plot ───────────────────────────────────────────────

def plot_expression_dotplot():
    path = OUTPUT_DIR / "celltype_expression.csv"
    if not path.exists():
        log.warning(f"  {path} not found — skipping expression dot plot")
        return
    expr_df = pd.read_csv(path)
    log.info("Plotting expression dot plot (fig1)...")

    mut_genes = [g for g in GWAS_LOCI if GWAS_LOCI[g]["subtype"] == "IDH-mut"]
    wt_genes  = [g for g in GWAS_LOCI if GWAS_LOCI[g]["subtype"] == "IDH-wt"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7),
                             gridspec_kw={"width_ratios": [2, 1]})
    norm = Normalize(vmin=0, vmax=expr_df["mean_expression"].quantile(0.95))

    for ax, genes, title, cmap_name in [
        (axes[0], mut_genes, "IDH-mutant risk loci",   "Reds"),
        (axes[1], wt_genes,  "IDH-wildtype risk loci", "Blues"),
    ]:
        cmap = matplotlib.colormaps[cmap_name]
        for i, gene in enumerate(genes):
            for j, ct in enumerate(EXPRESSION_CELL_TYPES):
                row = expr_df[(expr_df["gene"] == gene) &
                              (expr_df["cell_type"] == ct)]
                if len(row) == 0:
                    continue
                pct = row["pct_expressing"].iloc[0]
                me  = row["mean_expression"].iloc[0]
                if pd.notna(pct):
                    size = max(pct * 4, 5)
                elif pd.notna(me):
                    # HPA data: no pct_expressing, use mean_expression for size
                    size = max(me * 0.3, 5)
                else:
                    continue
                ax.scatter(j, len(genes) - 1 - i, s=size,
                           c=[cmap(norm(me))],
                           edgecolors="gray", linewidths=0.5, zorder=3)
        ax.set_xticks(range(len(EXPRESSION_CELL_TYPES)))
        ax.set_xticklabels(EXPRESSION_CELL_TYPES, rotation=45, ha="right",
                           fontsize=9)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(list(reversed(genes)), fontsize=10,
                           fontstyle="italic")
        ax.set_title(title, fontweight="bold")
        # Highlight oligodendrocyte lineage
        ax.axvspan(-0.5, 0.5, alpha=0.08, color="orange", zorder=1)
        ax.axvspan(0.5,  1.5, alpha=0.08, color="orange", zorder=1)
        ax.set_xlim(-0.5, len(EXPRESSION_CELL_TYPES) - 0.5)
        ax.grid(True, alpha=0.2)

    fig.suptitle("Expression of glioma GWAS genes across brain cell types"
                 "\n(orange shading = oligodendrocyte lineage)",
                 fontweight="bold")
    plt.tight_layout()
    out = OUTPUT_DIR / "fig1_expression_dotplot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved {out}")


# ── Fig 2: eQTL heatmap (BH-corrected) ──────────────────────────────────────

def plot_eqtl_heatmap():
    path = OUTPUT_DIR / "eqtl_corrected.csv"
    if not path.exists():
        log.warning(f"  {path} not found — skipping eQTL heatmap")
        return
    best = pd.read_csv(path)
    log.info("Plotting eQTL heatmap (fig2)...")

    cell_types = RANKING_CELL_TYPES        # no pseudobulk
    ct_labels  = ["OPC", "Oligo", "Astro", "ExN", "InN", "Micro"]

    # Locus ordering from summary
    summary_path = OUTPUT_DIR / "locus_summary.csv"
    if summary_path.exists():
        summary = pd.read_csv(summary_path)
        mut_order = summary[summary["subtype"] == "IDH-mut"].sort_values(
            "best_pvalue", na_position="last")["locus"].tolist()
        wt_order = summary[summary["subtype"] == "IDH-wt"].sort_values(
            "best_pvalue", na_position="last")["locus"].tolist()
    else:
        mut_order = [l for l in GWAS_LOCI if GWAS_LOCI[l]["subtype"] == "IDH-mut"]
        wt_order  = [l for l in GWAS_LOCI if GWAS_LOCI[l]["subtype"] == "IDH-wt"]
    locus_order = mut_order + wt_order

    M = np.full((len(locus_order), len(cell_types)), np.nan)
    P_raw = np.full_like(M, np.nan)
    P_bh  = np.full_like(M, np.nan)
    for i, locus in enumerate(locus_order):
        for j, ct in enumerate(cell_types):
            r = best[(best["locus"] == locus) & (best["cell_type"] == ct)]
            if len(r) == 0:
                continue
            p = r["pvalue"].iloc[0]
            b = r["beta"].iloc[0]
            M[i, j] = -np.log10(max(p, 1e-12)) * np.sign(b)
            P_raw[i, j] = p
            if "pvalue_bh" in r.columns:
                P_bh[i, j] = r["pvalue_bh"].iloc[0]

    fig, ax = plt.subplots(figsize=(9, 7))
    finite = M[np.isfinite(M)]
    vmax = max(np.abs(finite).max() if len(finite) else 3, 3)
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="nearest")

    ax.set_xticks(range(len(ct_labels)))
    ax.set_xticklabels(ct_labels, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(locus_order)))
    ax.set_yticklabels(locus_order, fontsize=10, fontstyle="italic")

    for tick, locus in zip(ax.get_yticklabels(), locus_order):
        tick.set_color(_subtype_color(GWAS_LOCI[locus]["subtype"]))

    # Annotate: BH-corrected significance
    for i in range(len(locus_order)):
        for j in range(len(cell_types)):
            p_bh = P_bh[i, j]
            p_raw = P_raw[i, j]
            if np.isnan(p_raw):
                ax.text(j, i, "-", ha="center", va="center",
                        fontsize=9, color="lightgray")
            elif not np.isnan(p_bh) and p_bh < 0.001:
                ax.text(j, i, "***", ha="center", va="center",
                        fontsize=10, fontweight="bold")
            elif not np.isnan(p_bh) and p_bh < 0.01:
                ax.text(j, i, "**", ha="center", va="center", fontsize=10)
            elif not np.isnan(p_bh) and p_bh < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=10)

    ax.axvspan(-0.5, 1.5, alpha=0.07, color="orange", zorder=0)
    ax.axhline(len(mut_order) - 0.5, color="black", linestyle="--",
               linewidth=1.5)

    plt.colorbar(im, ax=ax, label="-log10(p) × sign(β)", shrink=0.7)
    ax.set_title("Cell-type-specific eQTL effects (Bryois 2022)"
                 "\nStars = BH-corrected within locus  "
                 "(* <0.05  ** <0.01  *** <0.001)",
                 fontweight="bold", fontsize=10)
    plt.tight_layout()
    out = OUTPUT_DIR / "fig2_eqtl_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved {out}")


# ── Fig 3: Locus summary bar chart ──────────────────────────────────────────

def plot_locus_summary():
    path = OUTPUT_DIR / "locus_summary.csv"
    if not path.exists():
        log.warning(f"  {path} not found — skipping locus summary plot")
        return
    summary = pd.read_csv(path)
    log.info("Plotting locus summary (fig3)...")

    df = summary.dropna(subset=["best_pvalue"]).copy()
    df["neg_log_p"] = -np.log10(df["best_pvalue"])
    df = df.sort_values(["subtype", "neg_log_p"], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(range(len(df)), df["neg_log_p"],
                   color=[CT_COLORS.get(c, "#cccccc")
                          for c in df["best_cell_type"]],
                   edgecolor="black", linewidth=0.5)

    # Mark loci that don't survive BH correction
    if "best_survives_bh" in df.columns:
        for idx_bar, (_, row) in zip(range(len(df)), df.iterrows()):
            if not row.get("best_survives_bh", True):
                bars[idx_bar].set_hatch("//")
                bars[idx_bar].set_edgecolor("gray")

    ax.set_yticks(range(len(df)))
    labels = [f"{l} ({s})" for l, s in zip(df["locus"], df["subtype"])]
    ax.set_yticklabels(labels, fontsize=9)
    for tick, sub in zip(ax.get_yticklabels(), df["subtype"]):
        tick.set_color(_subtype_color(sub))

    ax.axvline(-np.log10(0.05), color="gray", linestyle="--", label="p=0.05")
    ax.axvline(-np.log10(0.001), color="gray", linestyle=":", label="p=0.001")
    ax.set_xlabel("-log10(strongest eQTL p-value)")
    ax.set_title("Strongest eQTL per locus (color = cell type)\n"
                 "Hatching = does not survive within-locus BH correction",
                 fontweight="bold", fontsize=10)

    handles = [mpatches.Patch(color=c, label=ct) for ct, c in CT_COLORS.items()]
    handles += [
        plt.Line2D([0], [0], color="gray", linestyle="--", label="p=0.05"),
        plt.Line2D([0], [0], color="gray", linestyle=":",  label="p=0.001"),
        mpatches.Patch(facecolor="white", edgecolor="gray",
                       hatch="//", label="Fails BH"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7)
    plt.tight_layout()
    out = OUTPUT_DIR / "fig3_locus_summary.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved {out}")


# ── Fig 4: Colocalization heatmap ────────────────────────────────────────────

def plot_coloc_heatmap():
    path = OUTPUT_DIR / "coloc_results.csv"
    if not path.exists():
        log.info("  No coloc results — skipping fig4")
        return
    coloc_df = pd.read_csv(path)
    if coloc_df.empty or "PP.H4" not in coloc_df.columns:
        log.info("  Empty coloc results — skipping fig4")
        return
    if coloc_df["PP.H4"].isna().all():
        log.info("  All PP.H4 are NaN — skipping fig4")
        return
    log.info("Plotting colocalization heatmap (fig4)...")

    cell_types = RANKING_CELL_TYPES
    ct_labels = ["OPC", "Oligo", "Astro", "ExN", "InN", "Micro"]

    locus_order = ([l for l in GWAS_LOCI if GWAS_LOCI[l]["subtype"] == "IDH-mut"]
                   + [l for l in GWAS_LOCI if GWAS_LOCI[l]["subtype"] == "IDH-wt"])

    # Show the BEST PP.H4 across all GWAS subtypes for each (locus, cell_type)
    M = np.full((len(locus_order), len(cell_types)), np.nan)
    for i, locus in enumerate(locus_order):
        for j, ct in enumerate(cell_types):
            r = coloc_df[(coloc_df["locus"] == locus) &
                         (coloc_df["cell_type"] == ct)]
            if len(r):
                valid = r["PP.H4"].dropna()
                if len(valid):
                    M[i, j] = valid.max()

    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(M, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=1, interpolation="nearest")

    ax.set_xticks(range(len(ct_labels)))
    ax.set_xticklabels(ct_labels, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(locus_order)))
    ax.set_yticklabels(locus_order, fontsize=10, fontstyle="italic")
    for tick, locus in zip(ax.get_yticklabels(), locus_order):
        tick.set_color(_subtype_color(GWAS_LOCI[locus]["subtype"]))

    # Annotate PP.H4 values
    for i in range(len(locus_order)):
        for j in range(len(cell_types)):
            v = M[i, j]
            if np.isnan(v):
                ax.text(j, i, "-", ha="center", va="center",
                        fontsize=8, color="lightgray")
            else:
                color = "white" if v > 0.6 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=8, color=color, fontweight="bold" if v >= COLOC_PP4 else "normal")

    # Separator
    n_mut = sum(1 for l in GWAS_LOCI if GWAS_LOCI[l]["subtype"] == "IDH-mut")
    ax.axhline(n_mut - 0.5, color="black", linestyle="--", linewidth=1.5)

    plt.colorbar(im, ax=ax, label="PP.H4 (shared causal variant)", shrink=0.7)
    ax.set_title("Colocalization PP.H4: GWAS × cell-type eQTL\n"
                 f"(bold ≥ {COLOC_PP4} = evidence for shared mechanism)",
                 fontweight="bold", fontsize=10)
    plt.tight_layout()
    out = OUTPUT_DIR / "fig4_coloc_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved {out}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    log.info("=" * 72)
    log.info("Generating figures")
    log.info("=" * 72)
    plot_expression_dotplot()
    plot_eqtl_heatmap()
    plot_locus_summary()
    plot_coloc_heatmap()
    log.info("Done.")


if __name__ == "__main__":
    main()
