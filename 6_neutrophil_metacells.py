import os
os.environ["NUMBA_DISABLE_JIT"] = "1"

import metacells as mc
import scanpy as sc
import numpy as np
import scipy.sparse as sp

base_in = "/home/jupyter/trem1/neutrophil_v2_012026/de_cellchat/trem1pos_neut/metacells_input"
base_out = "/home/jupyter/trem1/neutrophil_v2_012026/de_cellchat/trem1pos_neut/metacells_output"

os.makedirs(base_out, exist_ok=True)

seed = 123456
target_metacell_size = 50000

excluded_gene_names = [
    "IGHMBP2", "IGLL1", "IGLL5", "IGLON5", "NEAT1", "TMSB10", "TMSB4X"
]
excluded_gene_patterns = [
    "MT-.*", "ENS.*"
]

suspect_gene_names = [
    "PCNA", "MKI67", "TOP2A", "HIST1H1D", "FOS",
    "JUN", "HSP90AB1", "HSPA1A", "ISG15", "WARS"
]
suspect_gene_patterns = [
    "MCM[0-9]", "SMC[0-9]", "IFI.*"
]

def do_one(group_name, out_name):
    adata = sc.read_10x_mtx(os.path.join(base_in, group_name))

    mc.pl.analyze_clean_genes(
        adata,
        excluded_gene_names=excluded_gene_names,
        excluded_gene_patterns=excluded_gene_patterns,
        random_seed=seed,
    )
    mc.pl.pick_clean_genes(adata)

    mc.pl.analyze_clean_cells(
        adata,
        properly_sampled_min_cell_total=800,
        properly_sampled_max_cell_total=8000,
        properly_sampled_max_excluded_genes_fraction=0.1,
    )
    mc.pl.pick_clean_cells(adata)

    clean = mc.pl.extract_clean_data(adata)

    suspect_mask = mc.tl.find_named_genes(
        clean,
        names=suspect_gene_names,
        patterns=suspect_gene_patterns,
    )

    if sp.isspmatrix_csc(clean.X):
        clean.X = clean.X.tocsr()
    if "raw" in clean.layers and sp.isspmatrix_csc(clean.layers["raw"]):
        clean.layers["raw"] = clean.layers["raw"].tocsr()

    mc.pl.relate_genes(clean, random_seed=seed)

    modules = clean.var["related_genes_module"].to_numpy()
    suspect_modules = np.unique(modules[suspect_mask])
    suspect_modules = suspect_modules[suspect_modules >= 0]

    sim = mc.ut.get_vv_frame(clean, "related_genes_similarity")

    bad_modules = []
    for m in suspect_modules:
        m_mask = modules == m
        sim_block = sim.loc[m_mask, m_mask]
        if sim_block.mean().mean() > 0.75:
            bad_modules.append(m)

    forbidden_mask = suspect_mask.copy()
    for m in bad_modules:
        forbidden_mask |= (modules == m)

    forbidden_gene_names = sorted(clean.var_names[forbidden_mask])

    mc.pl.set_max_parallel_piles(1)

    mc.pl.divide_and_conquer_pipeline(
        clean,
        forbidden_gene_names=forbidden_gene_names,
        random_seed=seed,
        target_metacell_size=target_metacell_size,
    )

    out = clean.obs[["metacell"]].copy()
    out.index.name = "cell_barcode"
    out.to_csv(os.path.join(base_out, out_name))

print("Running Control")
do_one("Control", "control_metacells.csv")

print("Running TET2")
do_one("TET2", "tet2_metacells.csv")
