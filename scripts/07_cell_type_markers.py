"""
07_cell_type_markers.py
Identify cell types using marker genes
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white')

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"

# Load clustered data
print("Loading data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_clustered.h5ad")
print(f"Loaded: {adata.shape[0]} cells, {adata.obs['leiden'].nunique()} clusters")

# Define marker genes for immune cell types
markers = {
    'T cell (general)': ['CD3E', 'CD3D'],
    'CD8 T cell': ['CD8A', 'CD8B'],
    'CD4 T cell': ['CD4'],
    'T reg': ['FOXP3', 'IL2RA'],
    'NK cell': ['NKG7', 'GNLY', 'NCAM1'],
    'B cell': ['CD19', 'MS4A1', 'CD79A'],
    'Macrophage': ['CD68', 'CD14', 'FCGR3A'],
    'Dendritic cell': ['CLEC4C', 'IL3RA'],
    'Exhaustion': ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT'],
    'Cytotoxic': ['GZMB', 'PRF1', 'IFNG'],
}

# Flatten for dotplot
all_markers = []
for genes in markers.values():
    all_markers.extend(genes)

# Check which markers are in our data
markers_present = [g for g in all_markers if g in adata.var_names]
markers_missing = [g for g in all_markers if g not in adata.var_names]

print(f"\nMarkers found: {len(markers_present)}/{len(all_markers)}")
if markers_missing:
    print(f"Missing: {markers_missing}")

# Plot dotplot of markers by cluster
print("\nPlotting marker dotplot...")
sc.pl.dotplot(adata, markers_present, groupby='leiden', 
              standard_scale='var', show=False)
plt.savefig(f"{FIGURES_DIR}/06_marker_dotplot.png", dpi=150, bbox_inches='tight')
print(f"Saved: {FIGURES_DIR}/06_marker_dotplot.png")
plt.close()

# Plot markers on UMAP
print("\nPlotting markers on UMAP...")

# Key markers to visualize
key_markers = ['CD3E', 'CD8A', 'CD4', 'FOXP3', 'NKG7', 'CD68', 'CD19', 'PDCD1']
key_markers = [m for m in key_markers if m in adata.var_names]

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

for i, marker in enumerate(key_markers):
    sc.pl.umap(adata, color=marker, ax=axes[i], show=False, 
               title=marker, cmap='Reds')

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/07_marker_umap.png", dpi=150)
print(f"Saved: {FIGURES_DIR}/07_marker_umap.png")
plt.close()

# Print cluster sizes
print("\nCluster sizes:")
print(adata.obs['leiden'].value_counts().sort_index())

# Save
adata.write_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_clustered.h5ad")
print("\nDone!")
