#!/usr/bin/env python3
"""
01_laptop_query.py
------------------
Lightweight drop-in replacement for 01_expression_profiling.py.

Use this when the cluster lacks AVX2 (tiledbsoma crashes) or S3 access.
Run on any machine with internet and a modern CPU (laptop, workstation).

Queries ALL normal human brain cells from the Census (not pinned to a
specific dataset ID), making it robust to Census reorganisations.

Output: output/celltype_expression.csv  (same format as 01_expression_profiling.py)

Usage:
    pip install --only-binary=:all: cellxgene-census
    python3 01_laptop_query.py
    scp output/celltype_expression.csv you@cluster:.../pipeline/output/
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp

GENES = ['AKT3','C2orf80','CCDC26','CDKN2A','CDKN2B','EGFR','ETFA',
         'IDH1','LRIG1','MAML2','PHLDB1','RTEL1','TERT','ZBTB16']

CELL_TYPES = ["OPC","Oligodendrocyte","Astrocyte","Microglia",
              "Excitatory Neuron","Inhibitory Neuron","Endothelial","Pericyte/Vascular"]

CT_RULES = [
    (["oligodendrocyte precursor","opc"],        "OPC"),
    (["oligodendrocyte"],                        "Oligodendrocyte"),
    (["astro"],                                  "Astrocyte"),
    (["microglia"],                              "Microglia"),
    (["glutamatergic","excitatory"],             "Excitatory Neuron"),
    (["gabaergic","inhibitory","interneuron"],   "Inhibitory Neuron"),
    (["endothelial"],                            "Endothelial"),
    (["pericyte","vascular","smooth muscle"],    "Pericyte/Vascular"),
]

def map_ct(ct):
    if not isinstance(ct, str): return "Other"
    ctl = ct.lower()
    for kws, label in CT_RULES:
        if any(k in ctl for k in kws): return label
    return "Other"


def main():
    import cellxgene_census
    import tiledbsoma as soma

    os.environ["TILEDB_SOMA_INIT_BUFFER_BYTES"] = str(512 * 1024 * 1024)
    context = soma.SOMATileDBContext(tiledb_config={
        "vfs.s3.region": "us-west-2",
        "vfs.s3.no_sign_request": "true",
    })

    print("Opening CELLxGENE Census 2025-11-08...")
    census = cellxgene_census.open_soma(census_version="2025-11-08",
                                         context=context)

    print(f"Querying {len(GENES)} genes from all normal human brain tissue...")
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",
        var_value_filter=f"feature_name in {GENES}",
        obs_value_filter=(
            "tissue_general == 'brain' "
            "and is_primary_data == True "
            "and disease == 'normal'"
        ),
        obs_column_names=["cell_type", "tissue", "donor_id"],
    )
    census.close()
    print(f"Retrieved {adata.shape[0]:,} cells x {adata.shape[1]} genes")

    adata.obs["ct_major"] = adata.obs["cell_type"].apply(map_ct)
    n_other = (adata.obs["ct_major"] == "Other").sum()
    print(f"Cell-type mapping: {n_other:,}/{len(adata):,} mapped to Other")

    adata = adata[adata.obs["ct_major"].isin(CELL_TYPES)].copy()
    print(f"Kept {adata.shape[0]:,} cells in {len(CELL_TYPES)} major types")
    for ct, n in adata.obs["ct_major"].value_counts().items():
        print(f"  {ct:22s} {n:>9,}")

    gene_names = adata.var["feature_name"].tolist()
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)

    rows = []
    for i, gene in enumerate(gene_names):
        for ct in CELL_TYPES:
            mask = (adata.obs["ct_major"] == ct).values
            n = int(mask.sum())
            if n == 0:
                continue
            expr = X[mask, i]
            ne = int((expr > 0).sum())
            rows.append({
                "gene": gene,
                "cell_type": ct,
                "n_cells": n,
                "n_expressing": ne,
                "pct_expressing": round(100.0 * ne / n, 2),
                "mean_expression": round(float(expr.mean()), 4),
            })

    df = pd.DataFrame(rows)
    os.makedirs("output", exist_ok=True)
    df.to_csv("output/celltype_expression.csv", index=False)
    print(f"\nSaved output/celltype_expression.csv ({len(df)} rows)")

    pivot = df.pivot_table(index="gene", columns="cell_type",
                           values="pct_expressing")[CELL_TYPES].round(1)
    print(f"\n% expressing:\n{pivot.to_string()}")
    print("\nDone.")


if __name__ == "__main__":
    main()
