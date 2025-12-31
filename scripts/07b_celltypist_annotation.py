"""
07b_celltypist_annotation.py
Automated cell type annotation using CellTypist database
"""

import os
import scanpy as sc
import celltypist
from celltypist import models
import matplotlib.pyplot as plt

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"

# Load data
print("Loading data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_clustered.h5ad")
print(f"Loaded: {adata.shape[0]} cells")

# Download immune cell model (first time only)
print("\nDownloading immune cell model...")
models.download_models(model='Immune_All_Low.pkl')

# Load the model
model = models.Model.load(model='Immune_All_Low.pkl')
print(f"Model has {len(model.cell_types)} cell types")

# Run prediction
print("\nPredicting cell types (this may take a minute)...")
predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl', majority_voting=True)

# Add predictions to adata
adata = predictions.to_adata()

# Show results
print("\nPredicted cell types:")
print(adata.obs['majority_voting'].value_counts())

# Plot
print("\nPlotting...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='Clusters')
sc.pl.umap(adata, color='majority_voting', ax=axes[1], show=False, title='CellTypist annotation')

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/09_celltypist_annotation.png", dpi=150, bbox_inches='tight')
print(f"Saved: {FIGURES_DIR}/09_celltypist_annotation.png")
plt.close()

# Save
adata.write_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_celltypist.h5ad")
print(f"Saved: {PROCESSED_DIR}/GSE120575_melanoma_celltypist.h5ad")
