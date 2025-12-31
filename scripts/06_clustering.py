"""
06_clustering.py
Cluster cells and visualize
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white')

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"

# Load data with UMAP
print("Loading data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_pca_umap.h5ad")
print(f"Loaded: {adata.shape[0]} cells")

# Cluster at different resolutions
print("\nClustering...")
for res in [0.3, 0.5, 0.8]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    n = adata.obs[f'leiden_{res}'].nunique()
    print(f"  Resolution {res}: {n} clusters")

# Use 0.5 as default
adata.obs['leiden'] = adata.obs['leiden_0.5']

# Plot clusters
print("\nPlotting...")
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

sc.pl.umap(adata, color='leiden_0.3', ax=axes[0], show=False, 
           title=f'res=0.3 ({adata.obs["leiden_0.3"].nunique()} clusters)', legend_loc='on data')
sc.pl.umap(adata, color='leiden_0.5', ax=axes[1], show=False,
           title=f'res=0.5 ({adata.obs["leiden_0.5"].nunique()} clusters)', legend_loc='on data')
sc.pl.umap(adata, color='leiden_0.8', ax=axes[2], show=False,
           title=f'res=0.8 ({adata.obs["leiden_0.8"].nunique()} clusters)', legend_loc='on data')

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/05_clusters.png", dpi=150)
print(f"Saved: {FIGURES_DIR}/05_clusters.png")
plt.close()

# Save
adata.write_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_clustered.h5ad")
print(f"Saved: {PROCESSED_DIR}/GSE120575_melanoma_clustered.h5ad")
