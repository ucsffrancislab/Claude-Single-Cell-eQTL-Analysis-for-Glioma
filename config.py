#!/usr/bin/env python3
"""
Shared configuration for the glioma single-cell eQTL pipeline.

All locus definitions, cell-type mappings, file paths, and constants live here
so every script imports a single source of truth.
"""

import os
from pathlib import Path

# =============================================================================
# PATHS  (all relative to PIPELINE_DIR so the pipeline is relocatable)
# =============================================================================

PIPELINE_DIR = Path(os.environ.get("SCEQTL_PIPELINE_DIR",
                                    Path(__file__).resolve().parent))
DATA_DIR     = PIPELINE_DIR / "data"
BRYOIS_DIR   = DATA_DIR / "bryois"
GWAS_DIR     = DATA_DIR / "gwas"            # user-supplied GWAS summary stats
OUTPUT_DIR   = PIPELINE_DIR / "output"

for d in (DATA_DIR, BRYOIS_DIR, GWAS_DIR, OUTPUT_DIR):
    d.mkdir(parents=True, exist_ok=True)

# =============================================================================
# SIGNIFICANCE THRESHOLDS
# =============================================================================

NOMINAL_P    = 0.05          # per-test nominal threshold
BH_FDR       = 0.05          # Benjamini-Hochberg FDR for within-locus correction
COLOC_PP4    = 0.80          # posterior probability threshold for colocalization
N_PERM       = 10_000        # permutations for enrichment test

# =============================================================================
# CELLxGENE CENSUS
# =============================================================================

CENSUS_VERSION    = "2025-11-08"
ABC_ATLAS_DATASET = "6f7fd0f1-b445-4486-a377-6f07536f6c60"

# =============================================================================
# BRYOIS eQTL DATA
# =============================================================================

ZENODO_RECORD   = "7276971"
ZENODO_API      = f"https://zenodo.org/api/records/{ZENODO_RECORD}"
BRYOIS_N        = 196          # sample size in Bryois et al. 2022

# Map display name -> Zenodo file prefix
BRYOIS_CT_PREFIX = {
    "OPC":              "OPCs...COPs",
    "Astrocyte":        "Astrocytes",
    "ExcitatoryNeuron": "Excitatory.neurons",
    "InhibitoryNeuron": "Inhibitory.neurons",
    "Microglia":        "Microglia",
    "Oligodendrocyte":  "Oligodendrocytes",
}

# Pseudobulk kept for reference only — excluded from cell-type ranking
BRYOIS_PSEUDO_PREFIX = {"Pseudobulk": "pb"}

# All Bryois prefixes (for downloading)
BRYOIS_ALL_PREFIX = {**BRYOIS_CT_PREFIX, **BRYOIS_PSEUDO_PREFIX}

# Cell types that participate in "best cell type" ranking
# (Pseudobulk is a tissue-level aggregate, NOT a cell type)
RANKING_CELL_TYPES = list(BRYOIS_CT_PREFIX.keys())

# =============================================================================
# CELLxGENE EXPRESSION CELL TYPES  (superset, includes types not in Bryois)
# =============================================================================

EXPRESSION_CELL_TYPES = [
    "OPC", "Oligodendrocyte", "Astrocyte", "Microglia",
    "Excitatory Neuron", "Inhibitory Neuron",
    "Endothelial", "Pericyte/Vascular",
]

# =============================================================================
# GWAS LOCI  (Melin et al. 2017 Nat Genet, subtypes per Adel Fahmideh 2019)
# =============================================================================
#
# Each locus has:
#   rsid          – lead GWAS SNP
#   chr           – chromosome (int)
#   locus         – cytogenetic band
#   subtype       – IDH-mut or IDH-wt
#   pathway       – biological pathway
#   proxies       – [(rsid, r2), ...] when lead not in Bryois
#   alt_genes     – extra genes to test at this locus
#   gwas_beta     – log-OR from GWAS (OPTIONAL; set to None if unavailable)
#   gwas_se       – SE of log-OR      (OPTIONAL)
#   gwas_eaf      – effect allele freq (OPTIONAL)
#
# IMPORTANT: populate gwas_beta / gwas_se / gwas_eaf from your GWAS summary
# stats to enable risk-allele direction analysis and colocalization.

GWAS_LOCI = {
    # ---------- IDH-mutant associated ----------
    "CCDC26": {
        "rsid": "rs55705857", "chr": 8, "locus": "8q24.21",
        "subtype": "IDH-mut", "pathway": "MYC enhancer / lncRNA",
        # lead not in Bryois — use highest-r2 proxy only
        "proxies": [
            ("rs143893586", 0.636),
            ("rs55862293",  0.501),
            ("rs72714295",  0.494),
        ],
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "AKT3": {
        "rsid": "rs12076373", "chr": 1, "locus": "1q44",
        "subtype": "IDH-mut", "pathway": "PI3K/AKT",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "IDH1": {
        "rsid": "rs7572263", "chr": 2, "locus": "2q33.3",
        "subtype": "IDH-mut", "pathway": "Metabolic / G-CIMP",
        "alt_genes": ["C2orf80"],
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "LRIG1": {
        "rsid": "rs11706832", "chr": 3, "locus": "3p14.1",
        "subtype": "IDH-mut", "pathway": "EGFR negative regulator",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "MAML2": {
        "rsid": "rs7107785", "chr": 11, "locus": "11q21",
        "subtype": "IDH-mut", "pathway": "Notch co-activator",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "PHLDB1": {
        "rsid": "rs12803321", "chr": 11, "locus": "11q23.3",
        "subtype": "IDH-mut", "pathway": "PI3K/AKT activator",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "ETFA": {
        "rsid": "rs1801591", "chr": 15, "locus": "15q24.2",
        "subtype": "IDH-mut", "pathway": "Fatty acid oxidation",
        "proxies": [("rs77633900", None)],
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "ZBTB16": {
        "rsid": "rs648044", "chr": 11, "locus": "11q23.2",
        "subtype": "IDH-mut", "pathway": "Progenitor maintenance",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    # ---------- IDH-wildtype associated ----------
    "TERT": {
        "rsid": "rs10069690", "chr": 5, "locus": "5p15.33",
        "subtype": "IDH-wt", "pathway": "Telomere maintenance",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "EGFR": {
        "rsid": "rs75061358", "chr": 7, "locus": "7p11.2",
        "subtype": "IDH-wt", "pathway": "EGFR signaling",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "CDKN2A": {
        "rsid": "rs634537", "chr": 9, "locus": "9p21.3",
        "subtype": "IDH-wt", "pathway": "Cell cycle (p16/p14arf)",
        # lead not in Bryois; proxy r2 unknown — RESOLVE BEFORE PUBLICATION
        "proxies": [("rs613312", None)],
        "alt_genes": ["CDKN2B"],
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
    "RTEL1": {
        "rsid": "rs2297440", "chr": 20, "locus": "20q13.33",
        "subtype": "IDH-wt", "pathway": "Telomere maintenance / DNA repair",
        "gwas_beta": None, "gwas_se": None, "gwas_eaf": None,
    },
}

# Derived gene list
ALL_GENES = sorted(
    {g for g in GWAS_LOCI}
    | {ag for cfg in GWAS_LOCI.values() for ag in cfg.get("alt_genes", [])}
)

# Chromosomes we need
NEEDED_CHROMS = sorted({cfg["chr"] for cfg in GWAS_LOCI.values()})


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def best_proxy(cfg):
    """Return the single highest-r2 proxy for a locus, or None.

    If r2 is None (unknown), the proxy is kept but flagged.
    The lead SNP is always preferred when it has no proxies.
    """
    proxies = cfg.get("proxies", [])
    if not proxies:
        return None
    # Sort descending by r2; None sorts last
    ranked = sorted(proxies, key=lambda x: x[1] if x[1] is not None else -1,
                    reverse=True)
    return ranked[0]


def snps_for_locus(locus_name):
    """Return list of (snp_id, role, r2) to test for a locus.

    Strategy (fixes review item #2):
      - Always include the lead SNP.
      - If proxies exist, include ONLY the highest-r2 proxy
        (not all proxies — avoids inflating significance).
    """
    cfg = GWAS_LOCI[locus_name]
    snps = [(cfg["rsid"], "lead", 1.0)]
    bp = best_proxy(cfg)
    if bp is not None:
        snps.append((bp[0], "proxy", bp[1]))
    return snps


def genes_for_locus(locus_name):
    """Return list of gene symbols to test at a locus."""
    cfg = GWAS_LOCI[locus_name]
    return [locus_name] + cfg.get("alt_genes", [])


def setup_logging(name, log_file=None):
    """Consistent logging across pipeline scripts."""
    import sys, logging
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode="w"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
    )
    return logging.getLogger(name)
