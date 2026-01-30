import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

# --- 1. Configuration ---
CONFIG = {
    'RNA_DIR': '/s1/chengd/data/mix_hg38_mm10/outs/filtered_feature_bc_matrix/',
    'OUT_DIR': '/s1/chengd/project/multiome/02_cellline/species_separation/',
    'QC': {
        'min_genes': 300, 
        'min_counts': 500, 
        'max_mito': 30
    }
}
os.makedirs(CONFIG['OUT_DIR'], exist_ok=True)

# --- 2. Core Functions ---

def load_rna_data(path):
    """Load RNA data and perform initial Quality Control"""
    print(f"\n[Loading RNA] from: {path}")
    
    # Load 10x MTX format
    adata = sc.read_10x_mtx(path, var_names='gene_ids', cache=False)
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Calculate Mitochondrial percentage
    mito_genes = adata.var_names.str.contains('^MT-|^mt-', regex=True)
    if mito_genes.sum() > 0:
        # Calculate ratio of mito counts to total counts
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
        mito_counts = np.array(adata[:, mito_genes].X.sum(axis=1)).flatten()
        adata.obs['percent_mito'] = (mito_counts / (total_counts + 1e-10)) * 100
    else:
        adata.obs['percent_mito'] = 0
    
    # Apply filtering thresholds
    mask = (adata.obs.n_genes_by_counts > CONFIG['QC']['min_genes']) & \
           (adata.obs.total_counts > CONFIG['QC']['min_counts']) & \
           (adata.obs.percent_mito < CONFIG['QC']['max_mito'])
    
    adata_filtered = adata[mask].copy()
    print(f"   -> {adata_filtered.n_obs} cells remaining after QC")
    return adata_filtered

def classify_species(adata):
    """Classify species using KDE-based automatic thresholding for RNA"""
    # Identify Human vs Mouse genes based on prefix
    feats = adata.var_names
    is_hs = feats.str.startswith(('GRCh38_', 'hg19_'))
    is_mm = feats.str.startswith(('mm10_', 'mm9_'))

    # Calculate counts per species
    X = adata.X
    count_hs = np.array(X[:, is_hs].sum(axis=1)).flatten()
    count_mm = np.array(X[:, is_mm].sum(axis=1)).flatten()
    ratio_hs = count_hs / (count_hs + count_mm + 1e-10)

    # Thresholding logic: Use Kernel Density Estimation (KDE) to find the valley
    kde = stats.gaussian_kde(ratio_hs)
    x_range = np.linspace(0, 1, 1000)
    y_range = kde(x_range)
    
    # Search for valley in the middle range (0.15 to 0.85)
    mid_mask = (x_range > 0.15) & (x_range < 0.85)
    if mid_mask.sum() > 0:
        valley = x_range[mid_mask][np.argmin(y_range[mid_mask])]
    else:
        valley = 0.15 # Fallback
    
    low_thresh, high_thresh = valley, 1 - valley

    # Assign species labels
    species = np.full(len(adata), 'Doublet', dtype=object)
    species[ratio_hs > high_thresh] = 'Human'
    species[ratio_hs < low_thresh] = 'Mouse'
    
    adata.obs['count_hs'] = count_hs
    adata.obs['count_mm'] = count_mm
    adata.obs['ratio_hs'] = ratio_hs
    adata.obs['species'] = species

    # Summary Statistics
    stats_dict = {
        'n_human': np.sum(species == 'Human'),
        'n_mouse': np.sum(species == 'Mouse'),
        'n_doublet': np.sum(species == 'Doublet'),
        'rate': np.sum(species == 'Doublet') / len(species) * 100
    }
    
    print(f"   [RNA] Results: Human={stats_dict['n_human']}, Mouse={stats_dict['n_mouse']}, Doublet={stats_dict['n_doublet']}")
    print(f"   [RNA] Doublet Rate: {stats_dict['rate']:.2f}%")
    return adata, stats_dict

def plot_species_mixing(adata, stats_res):
    """Generate Barnyard plot for species separation"""
    plt.figure(figsize=(8, 7), dpi=300)
    colors = {'Human': '#E74C3C', 'Mouse': '#377EB8', 'Doublet': '#984EA3'}
    
    df = adata.obs
    for sp in ['Doublet', 'Mouse', 'Human']:
        mask = df['species'] == sp
        plt.scatter(
            df.loc[mask, 'count_hs'], 
            df.loc[mask, 'count_mm'], 
            c=colors[sp], 
            s=15, 
            alpha=0.6, 
            edgecolors='none',
            label=f"{sp}: {mask.sum()} ({mask.sum()/len(df)*100:.1f}%)"
        )
    
    # Formatting the plot
    limit = max(df['count_hs'].max(), df['count_mm'].max()) * 1.05
    plt.plot([0, limit], [0, limit], 'k--', alpha=0.2)
    plt.title(f"RNA-seq Species Mixing\n(Doublet Rate: {stats_res['rate']:.2f}%)", fontsize=14, fontweight='bold')
    plt.xlabel("Human UMIs")
    plt.ylabel("Mouse UMIs")
    plt.legend(frameon=True, fontsize=10)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    
    plt.tight_layout()
    save_path = os.path.join(CONFIG['OUT_DIR'], 'RNA_Species_Mixing_Barnyard.pdf')
    plt.savefig(save_path)
    print(f"\n✅ Plot saved to: {save_path}")

# --- Main Execution ---

if __name__ == "__main__":
    # 1. Load and Filter
    adata_rna = load_rna_data(CONFIG['RNA_DIR'])
    
    # 2. Classify Species
    adata_rna, stats_rna = classify_species(adata_rna)
    
    # 3. Visualization
    plot_species_mixing(adata_rna, stats_rna)
