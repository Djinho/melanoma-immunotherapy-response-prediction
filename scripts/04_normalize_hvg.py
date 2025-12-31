"""
04_normalize_hvg.py
Normalize data and find highly variable genes
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.verbosity = 2

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"

# Load filtered data
print("Loading filtered data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_filtered.h5ad")
print(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Store raw counts
print("\nStoring raw counts in layer...")
adata.layers['counts'] = adata.X.copy()

# Normalize
print("Normalizing to 10,000 counts per cell...")
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform
print("Log transforming...")
sc.pp.log1p(adata)

# Store normalized
adata.layers['normalized'] = adata.X.copy()

# Find highly variable genes
print("\nFinding highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3', layer='counts')
n_hvg = adata.var['highly_variable'].sum()
print(f"Found {n_hvg} highly variable genes")

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
# Plot
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(f"{FIGURES_DIR}/02_highly_variable_genes.png", dpi=150, bbox_inches='tight')
print(f"Saved: {FIGURES_DIR}/02_highly_variable_genes.png")
plt.close()
# Save
output = f"{PROCESSED_DIR}/GSE120575_melanoma_normalized.h5ad"
adata.write_h5ad(output)
print(f"\nSaved: {output}")
