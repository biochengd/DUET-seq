#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DUET-SEQ ANALYSIS PIPELINE
Script: 03_rna_velocity_analysis_scvelo.py
Description: RNA velocity analysis of Spermatogonia (SPG) using scVelo dynamical mode.
             Includes quality filtering based on likelihood and latent time correction.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

# Global Plotting Settings
scv.set_figure_params(style='white', frameon=False, color_map='gnuplot', vector='wheel')
plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['savefig.dpi'] = 300

# 1. Configuration and Data Loading --------------------------------------------

def prepare_ann_data(data_dir):
    """Load matrices exported from R and construct AnnData object."""
    # Load matrices
    spliced = pd.read_csv(f'{data_dir}/spliced_matrix.csv', index_col=0)
    unspliced = pd.read_csv(f'{data_dir}/unspliced_matrix.csv', index_col=0)
    metadata = pd.read_csv(f'{data_dir}/metadata.csv', index_col=0)
    umap_coords = pd.read_csv(f'{data_dir}/umap_coordinates.csv', index_col=0)
    
    # Create AnnData (transpose to Cell x Gene)
    adata = sc.AnnData(X=spliced.T)
    adata.layers['spliced'] = spliced.T.values
    adata.layers['unspliced'] = unspliced.T.values
    
    # Add metadata and embeddings
    adata.obs = metadata.loc[adata.obs_names].copy()
    adata.obsm['X_umap'] = umap_coords.loc[adata.obs_names].values
    
    # Clean category columns
    if 'celltype_spg_lvl2' in adata.obs.columns:
        adata.obs['celltype'] = adata.obs['celltype_spg_lvl2'].astype('category')
        
    return adata

# Define Anonymized Input Path
INPUT_DIR = "./data/velocity_inputs"
OUTPUT_DIR = "./output/velocity"
os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = prepare_ann_data(INPUT_DIR)

# 2. Preprocessing Pipeline ----------------------------------------------------

# Basic filtering and normalization
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# Compute moments (required for velocity)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 3. Dynamical Velocity Model --------------------------------------------------

# Recover dynamics (fits the kinetic splicing model)
scv.tl.recover_dynamics(adata, n_jobs=16)

# Calculate velocity and velocity graph
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata, n_jobs=16)

# 4. Advanced Quality Filtering (Likelihood Filtering) -------------------------

def filter_by_likelihood(adata, threshold=0.05):
    """Retain only genes with high likelihood scores from the dynamical model."""
    mask = adata.var['fit_likelihood'] > threshold
    genes_to_keep = adata.var_names[mask].tolist()
    adata_filtered = adata[:, genes_to_keep].copy()
    
    # Recompute graph on high-quality subset
    scv.tl.velocity_graph(adata_filtered, n_jobs=16)
    return adata_filtered

# Apply high-quality gene filter
adata_hq = filter_by_likelihood(adata, threshold=0.02)

# 5. Latent Time and Directionality Correction ---------------------------------

# Define the biological root (SSC) to anchor velocity direction
if 'SPG_SSC' in adata_hq.obs['celltype'].values:
    # Identify the first SSC cell as the root
    root_idx = np.where(adata_hq.obs['celltype'] == 'SPG_SSC')[0][0]
    adata_hq.uns['iroot'] = root_idx
    
    # Compute Latent Time based on root cell
    scv.tl.latent_time(adata_hq)
    
    # Refine graph based on latent time
    scv.tl.velocity_graph(adata_hq, n_jobs=16)

# 6. Visualization -------------------------------------------------------------

# Stream embedding plot
scv.pl.velocity_embedding_stream(
    adata_hq, 
    basis="X_umap", 
    color="celltype", 
    palette=sns.color_palette("husl", len(adata_hq.obs['celltype'].unique())),
    title="SPG RNA Velocity Stream (Dynamical Model)",
    save=f"{OUTPUT_DIR}/spg_velocity_stream.pdf"
)

# Latent Time plot (Heatmap showing progression from SSC)
scv.pl.velocity_embedding(
    adata_hq, 
    basis="X_umap", 
    arrow_length=0, 
    color="latent_time", 
    color_map="gnuplot",
    title="SPG Latent Time",
    save=f"{OUTPUT_DIR}/spg_latent_time.pdf"
)

# 7. Driver Gene Identification ------------------------------------------------

# Rank genes by velocity correlation within specific cell types
scv.tl.rank_velocity_genes(adata_hq, groupby='celltype', min_corr=.3)
velocity_genes_df = pd.DataFrame(adata_hq.uns['rank_velocity_genes']['names'])
velocity_genes_df.to_csv(f"{OUTPUT_DIR}/ranked_velocity_genes.csv")

# Plot top driver genes for SSC -> Differentiation
top_genes = velocity_genes_df['SPG_SSC'].head(5).tolist()
scv.pl.velocity(
    adata_hq, 
    top_genes, 
    ncols=2, 
    basis="X_umap", 
    color="celltype",
    save=f"{OUTPUT_DIR}/top_ssc_driver_genes.png"
)

# 8. Save Final Object ---------------------------------------------------------
adata_hq.write(f"{OUTPUT_DIR}/spg_velocity_processed.h5ad")
