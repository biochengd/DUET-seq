#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DUET-SEQ ANALYSIS PIPELINE
Script: 04_celloracle_grn_perturbation_analysis.py
Description: Construction of base GRNs using Cicero and TF motifs, 
             followed by in silico perturbation simulations (KO) 
             for stage-specific trajectory analysis.
"""

import os
import gc
import time
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.applications import Gradient_calculator, Oracle_development_module, Oracle_systematic_analysis_helper

# 1. Configuration and Paths ---------------------------------------------------

# Set working directory to local data folder
DATA_DIR = "./data/celloracle"
OUTS_DIR = "./output/celloracle"
GENOME_DIR = "./data/genome/mm10"
REF_GENOME = "mm10"

os.makedirs(OUTS_DIR, exist_ok=True)

# Parameters
COACCESS_THRESHOLD = 0.8
N_GRID = 40
MIN_MASS = 18
ALPHA = 10
N_NEIGHBORS = 200

# 2. Base GRN Construction (Cicero + Motifs) ----------------------------------

def construct_base_grn():
    """Build the base Gene Regulatory Network using Cicero and Motif scanning."""
    # Load peaks and Cicero connections
    all_peaks_df = pd.read_csv(f"{DATA_DIR}/all_peaks_nc_spg.csv")
    peaks = np.array([p.replace('-', '_') for p in all_peaks_df.x.values])
    
    cicero_conn = pd.read_csv(f"{DATA_DIR}/cicero_connections_nc_spg.csv")
    cicero_conn['Peak1'] = cicero_conn['Peak1'].str.replace('-', '_')
    cicero_conn['Peak2'] = cicero_conn['Peak2'].str.replace('-', '_')

    # Annotate TSS and integrate Cicero
    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=REF_GENOME)
    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, 
                                                   cicero_connections=cicero_conn)

    # Filter by co-accessibility
    integrated_filtered = integrated[integrated.coaccess >= COACCESS_THRESHOLD].copy()

    # Initialize TFinfo and perform Motif Scan (or load pre-scanned data)
    # Note: Scanning is computationally intensive; usually loaded from file.
    tfi = ma.load_TFinfo(file_path=f"{DATA_DIR}/spg.celloracle.tfinfo")
    
    # Filter motifs by score and create base GRN
    tfi.filter_motifs_by_score(threshold=10)
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    base_grn = tfi.to_dataframe()
    base_grn.to_parquet(f"{DATA_DIR}/base_GRN_dataframe_spg.parquet")
    return base_grn

# 3. Oracle Object Initialization ----------------------------------------------

def setup_oracle_object(adata_path, base_grn):
    """Instantiate and prepare the Oracle object with scRNA-seq data."""
    adata = sc.read(adata_path)
    
    # Ensure TF list from motif analysis matches genes in AnnData
    with open(f'{DATA_DIR}/tf_direct_spg.txt', 'r') as f:
        tf_list = [line.strip() for line in f]
    
    # Filter AnnData for HVGs + TFs
    hvgs = adata.var.index[adata.var['highly_variable']].tolist()
    all_genes = sorted(list(set(hvgs).union(set(tf_list))))
    all_genes = [g for g in all_genes if g in adata.var_names]
    adata = adata[:, all_genes].copy()
    
    sc.pp.normalize_per_cell(adata)

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, 
                                       cluster_column_name="celltype_spg_lvl2", 
                                       embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=base_grn)
    
    # KNN Imputation
    oracle.perform_PCA()
    n_comps = 30 # Selected based on variance elbow
    k = int(0.025 * adata.shape[0])
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, n_jobs=4)
    
    return oracle

# 4. GRN Fitting and Single Gene Simulation ------------------------------------

def run_single_gene_perturbation(oracle, gene_of_interest="Crem"):
    """Simulate KO effect for a single gene and calculate vector fields."""
    # Fit GRN for simulation
    links = oracle.get_links(cluster_name_for_GRN_unit="celltype_spg_lvl2", alpha=ALPHA)
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=ALPHA, use_cluster_specific_TFdict=True)

    # Simulate KO (Set expression to 0)
    oracle.simulate_shift(perturb_condition={gene_of_interest: 0.0}, n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=N_NEIGHBORS, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # Calculate Gradient/Pseudotime Flow
    gradient = Gradient_calculator(oracle_object=oracle, 
                                   pseudotime_key="palantir_pseudotime")
    gradient.calculate_p_mass(smooth=0.8, n_grid=N_GRID, n_neighbors=N_NEIGHBORS)
    gradient.calculate_mass_filter(min_mass=MIN_MASS)
    gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly": 7})
    gradient.calculate_gradient()
    
    return oracle, gradient, links

# 5. Systematic/Batch KO Analysis ----------------------------------------------

def run_batch_simulation(oracle_path, links_path, gradient_obj, tf_list_path):
    """Perform KO simulations for a list of core TFs (e.g., FigR-identified)."""
    
    # Load core TFs (e.g., intersection of FigR and CellOracle)
    core_tfs = pd.read_csv(tf_list_path).iloc[:, 0].tolist()
    
    output_hdf5 = f"{DATA_DIR}/batch_ko_results.celloracle.hdf5"
    index_dict = {
        "Whole_cells": None,
        "SPG_SSC": np.where(sc.read(oracle_path).obs["celltype"] == "SPG_SSC")[0]
        # Add other celltype indices as needed
    }

    def reload_fresh_oracle():
        """Memory optimization: Reload oracle to prevent memory leak across simulations."""
        obj = co.load_hdf5(oracle_path)
        lnks = co.load_links(links_path)
        lnks.filter_links()
        obj.get_cluster_specific_TFdict_from_Links(links_object=lnks)
        obj.fit_GRN_for_simulation(alpha=ALPHA, use_cluster_specific_TFdict=True)
        return obj

    for tf in core_tfs:
        print(f"Simulating KO for: {tf}")
        try:
            current_oracle = reload_fresh_oracle()
            if tf not in current_oracle.active_regulatory_genes:
                continue
            
            # Simulate
            current_oracle.simulate_shift(perturb_condition={tf: 0}, n_propagation=3)
            current_oracle.estimate_transition_prob(n_neighbors=N_NEIGHBORS)
            current_oracle.calculate_embedding_shift(sigma_corr=0.05)

            # Store results in development module
            for name, idx in index_dict.items():
                dev = Oracle_development_module()
                dev.load_differentiation_reference_data(gradient_object=gradient_obj)
                dev.load_perturb_simulation_data(oracle_object=current_oracle, cell_idx_use=idx, name=name)
                dev.calculate_inner_product()
                dev.calculate_digitized_ip(n_bins=10)
                dev.set_hdf_path(path=output_hdf5) 
                dev.dump_hdf5(gene=tf, misc=name)
            
            del current_oracle
            gc.collect()
        except Exception as e:
            print(f"Error simulating {tf}: {e}")

# 6. Result Visualization ------------------------------------------------------

def plot_perturbation_scores(hdf5_path):
    """Compare perturbation scores (PS) between stages (SSC vs Meiotic)."""
    helper = Oracle_systematic_analysis_helper(hdf5_file_path=hdf5_path)
    ps_me = helper.calculate_negative_ps_p_value(misc="SPG_Meiotic")
    ps_ssc = helper.calculate_negative_ps_p_value(misc="SPG_SSC")

    # Log transform and merge
    ps_me['ps_sum'] = np.log1p(ps_me['ps_sum'])
    ps_ssc['ps_sum'] = np.log1p(ps_ssc['ps_sum'])
    
    ps_merged = pd.merge(ps_me, ps_ssc, on="gene", suffixes=('_meiotic', '_ssc'))
    
    # Scatter plot
    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=ps_merged, x='ps_sum_ssc', y='ps_sum_meiotic', color='grey', s=10)
    
    # Highlight specific drivers (Example: Dazl)
    highlight = ['Dazl']
    for gene in highlight:
        if gene in ps_merged['gene'].values:
            point = ps_merged[ps_merged['gene'] == gene]
            plt.scatter(point.ps_sum_ssc, point.ps_sum_meiotic, color='red', s=40)
            plt.text(point.ps_sum_ssc, point.ps_sum_meiotic, gene, color='red')

    plt.title("In Silico Perturbation Score Comparison")
    plt.xlabel("Perturbation Score [SPG_SSC]")
    plt.ylabel("Perturbation Score [SPG_Meiotic]")
    plt.savefig(f"{OUTS_DIR}/ps_comparison_ssc_vs_meiotic.pdf")
    plt.show()

# Main Execution Flow
if __name__ == "__main__":
    # 1. Build Base GRN
    base_grn_df = construct_base_grn()
    
    # 2. Setup Oracle
    oracle_obj = setup_oracle_object(f"{DATA_DIR}/spg_counts.h5ad", base_grn_df)
    
    # 3. Single Gene Test and Reference Gradient
    oracle_obj, grad_obj, links_obj = run_single_gene_perturbation(oracle_obj, "Crem")
    
    # 4. Systematic Analysis
    run_batch_simulation(f"{DATA_DIR}/spg.oracle", f"{DATA_DIR}/spg.links", grad_obj, f"{DATA_DIR}/core_tfs.csv")
    
    # 5. Summary Plots
    plot_perturbation_scores(f"{DATA_DIR}/batch_ko_results.celloracle.hdf5")
