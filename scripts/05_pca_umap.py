"""
05_pca_umap.py
Dimensionality reduction: PCA and UMAP
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white')

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"

# Load normalized data
print("Loading normalized data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_normalized.h5ad")
print(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"HVGs: {adata.var['highly_variable'].sum()}")

# Scale data (only HVGs)
print("\nScaling HVGs...")
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)

# PCA
print("Running PCA...")
sc.tl.pca(adata_hvg, n_comps=50)

# Copy results back
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']

# Plot variance explained
print("\nPlotting PCA variance...")
variance = adata.uns['pca']['variance_ratio']
cumvar = np.cumsum(variance)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(range(1, 51), variance, 'o-')
axes[0].set_xlabel('PC')
axes[0].set_ylabel('Variance explained')
axes[0].set_title('Scree plot')

axes[1].plot(range(1, 51), cumvar, 'o-')
axes[1].axhline(0.9, color='red', linestyle='--', label='90%')
axes[1].set_xlabel('PC')
axes[1].set_ylabel('Cumulative variance')
axes[1].legend()

n_pcs_90 = np.argmax(cumvar >= 0.9) + 1
print(f"PCs needed for 90% variance: {n_pcs_90}")

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/03_pca_variance.png", dpi=150)
print(f"Saved: {FIGURES_DIR}/03_pca_variance.png")
plt.close()

# Compute neighbors
print("\nComputing neighbors...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# UMAP
print("Running UMAP...")
sc.tl.umap(adata)

# Plot UMAP colored by metadata
print("\nPlotting UMAP...")
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

sc.pl.umap(adata, color='response', ax=axes[0], show=False, title='Response')
sc.pl.umap(adata, color='timepoint', ax=axes[1], show=False, title='Timepoint')
sc.pl.umap(adata, color='therapy', ax=axes[2], show=False, title='Therapy')

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/04_umap_metadata.png", dpi=150)
print(f"Saved: {FIGURES_DIR}/04_umap_metadata.png")
plt.close()

# Save
output = f"{PROCESSED_DIR}/GSE120575_melanoma_pca_umap.h5ad"
adata.write_h5ad(output)
print(f"\nSaved: {output}")
