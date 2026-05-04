#!/usr/bin/env python3
"""
01_expression_profiling.py
--------------------------
Query CZ CELLxGENE Census for per-cell-type expression of all glioma GWAS genes.

Requires network access to public S3 (us-west-2).
Run on a node with internet — login node is fine.

Output: output/celltype_expression.csv
"""

import sys, os
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Filter only the specific noisy warnings from tiledb/census
warnings.filterwarnings("ignore", category=FutureWarning, module="tiledbsoma")
warnings.filterwarnings("ignore", category=FutureWarning, module="cellxgene_census")
warnings.filterwarnings("ignore", message=".*X_name.*")

from config import (
    CENSUS_VERSION, ABC_ATLAS_DATASET, ALL_GENES,
    EXPRESSION_CELL_TYPES, OUTPUT_DIR, setup_logging,
)

log = setup_logging("01-expr", OUTPUT_DIR / "01_expression.log")


# ── Cell-type mapping ────────────────────────────────────────────────────────
# CELLxGENE uses fine-grained ontology labels.  Collapse to 8 major types.
# Ordering matters: check "oligodendrocyte precursor" before "oligodendrocyte".

_CT_RULES = [
    (["oligodendrocyte precursor", "opc"],              "OPC"),
    (["oligodendrocyte"],                               "Oligodendrocyte"),
    (["astro"],                                         "Astrocyte"),
    (["microglia"],                                     "Microglia"),
    (["glutamatergic", "excitatory"],                   "Excitatory Neuron"),
    (["gabaergic", "inhibitory", "interneuron"],        "Inhibitory Neuron"),
    (["endothelial"],                                   "Endothelial"),
    (["pericyte", "vascular", "smooth muscle"],         "Pericyte/Vascular"),
]

def map_cell_type(ct):
    if not isinstance(ct, str):
        return "Other"
    ctl = ct.lower()
    for keywords, label in _CT_RULES:
        if any(kw in ctl for kw in keywords):
            return label
    return "Other"


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    import cellxgene_census
    import tiledbsoma as soma

    os.environ["TILEDB_SOMA_INIT_BUFFER_BYTES"] = str(512 * 1024 * 1024)

    context = soma.SOMATileDBContext(tiledb_config={
        "vfs.s3.region": "us-west-2",
        "vfs.s3.no_sign_request": "true",
    })

    log.info(f"Opening CELLxGENE Census {CENSUS_VERSION}")
    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION,
                                         context=context)

    log.info(f"Querying {len(ALL_GENES)} genes: {', '.join(ALL_GENES)}")
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",
        var_value_filter=f"feature_name in {ALL_GENES}",
        obs_value_filter=(
            f"dataset_id == '{ABC_ATLAS_DATASET}' "
            f"and is_primary_data == True "
            f"and disease == 'normal'"
        ),
        obs_column_names=["cell_type", "assay", "tissue", "donor_id"],
    )
    census.close()
    log.info(f"Retrieved {adata.shape[0]:,} cells x {adata.shape[1]} genes")

    # Map to major cell types
    adata.obs["cell_type_major"] = adata.obs["cell_type"].apply(map_cell_type)
    n_other = (adata.obs["cell_type_major"] == "Other").sum()
    n_total = len(adata)
    log.info(f"Cell-type mapping: {n_other:,}/{n_total:,} "
             f"({100*n_other/n_total:.1f}%) mapped to 'Other'")
    assert n_other / n_total < 0.10, (
        f">{10}% cells unmapped — check CELLxGENE ontology labels"
    )

    mask = adata.obs["cell_type_major"].isin(EXPRESSION_CELL_TYPES)
    adata = adata[mask].copy()
    log.info(f"Kept {adata.shape[0]:,} cells in {len(EXPRESSION_CELL_TYPES)} types")
    for ct, n in adata.obs["cell_type_major"].value_counts().items():
        log.info(f"  {ct:22s} {n:>9,}")

    # Compute per-(gene, cell_type) statistics
    gene_names = adata.var["feature_name"].tolist()
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)

    rows = []
    for i, gene in enumerate(gene_names):
        for ct in EXPRESSION_CELL_TYPES:
            ct_mask = (adata.obs["cell_type_major"] == ct).values
            n_cells = int(ct_mask.sum())
            if n_cells == 0:
                continue
            expr = X[ct_mask, i]
            n_expr = int((expr > 0).sum())
            rows.append({
                "gene":            gene,
                "cell_type":       ct,
                "n_cells":         n_cells,
                "n_expressing":    n_expr,
                "pct_expressing":  round(100.0 * n_expr / n_cells, 2),
                "mean_expression": round(float(expr.mean()), 4),
            })

    expr_df = pd.DataFrame(rows)
    out = OUTPUT_DIR / "celltype_expression.csv"
    expr_df.to_csv(out, index=False)
    log.info(f"Saved {out} ({len(expr_df)} rows)")

    # Pretty-print pivot
    pivot = expr_df.pivot_table(
        index="gene", columns="cell_type", values="pct_expressing"
    )[EXPRESSION_CELL_TYPES].round(1)
    log.info(f"\n% expressing:\n{pivot.to_string()}\n")
    log.info("Done.")


if __name__ == "__main__":
    main()
