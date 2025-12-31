"""
03_qc_filter.py
Quality control and filtering only
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

# Setup
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white')

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# Load data
print("Loading data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma.h5ad")
print(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Calculate QC metrics
print("\nCalculating QC metrics...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

print(f"Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"Median MT%: {adata.obs['pct_counts_mt'].median():.1f}%")

# Plot QC
print("\nPlotting QC metrics...")
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

axes[0].hist(adata.obs['n_genes_by_counts'], bins=50)
axes[0].axvline(200, color='red', linestyle='--')
axes[0].set_xlabel('Genes per cell')
axes[0].set_title('Filter: keep > 200')

axes[1].hist(adata.obs['pct_counts_mt'], bins=50)
axes[1].axvline(20, color='red', linestyle='--')
axes[1].set_xlabel('MT %')
axes[1].set_title('Filter: keep < 20%')

axes[2].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], s=1, alpha=0.3)
axes[2].set_xlabel('Total counts')
axes[2].set_ylabel('Genes detected')

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/01_qc_metrics.png", dpi=150)
print(f"Saved: {FIGURES_DIR}/01_qc_metrics.png")
plt.close()

# Filter
print("\nFiltering...")
print(f"Before: {adata.shape[0]} cells")

sc.pp.filter_cells(adata, min_genes=200)
print(f"After min_genes=200: {adata.shape[0]} cells")

adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()
print(f"After MT < 20%: {adata.shape[0]} cells")

sc.pp.filter_genes(adata, min_cells=10)
print(f"After min_cells=10: {adata.shape[1]} genes")

# Save
output = f"{PROCESSED_DIR}/GSE120575_melanoma_filtered.h5ad"
adata.write_h5ad(output)
print(f"\nSaved: {output}")
print(f"Final: {adata.shape[0]} cells × {adata.shape[1]} genes")
