"""
09_liana_communication.py
Cell-cell communication analysis using LIANA
"""

import scanpy as sc
import liana as li
import matplotlib.pyplot as plt
import pandas as pd
import warnings

warnings.filterwarnings('ignore')

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

# Load pre-treatment data
print("Loading pre-treatment data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_pre_treatment.h5ad")
print(f"Cells: {adata.shape[0]}")

# Setup
adata.obs['cell_type'] = adata.obs['majority_voting']

# Filter cell types with ≥30 cells
cell_counts = adata.obs['cell_type'].value_counts()
valid_types = cell_counts[cell_counts >= 30].index.tolist()
print(f"Cell types with ≥30 cells: {len(valid_types)}")

adata_filtered = adata[adata.obs['cell_type'].isin(valid_types)].copy()
print(f"Cells after filtering: {adata_filtered.shape[0]}")

# Run LIANA
print("\nRunning LIANA...")
li.mt.rank_aggregate(
    adata_filtered,
    groupby='cell_type',
    resource_name='consensus',
    use_raw=False,
    verbose=True
)

# Get results
liana_results = adata_filtered.uns['liana_res']
print(f"\nTotal interactions: {len(liana_results)}")

# Save full results
liana_results.to_csv(f"{RESULTS_DIR}/liana_all_interactions.csv", index=False)
print(f"Saved: {RESULTS_DIR}/liana_all_interactions.csv")

# Show top 20
print("\nTop 20 interactions:")
top = liana_results.sort_values('magnitude_rank').head(20)
print(top[['source', 'target', 'ligand_complex', 'receptor_complex']].to_string())

# Manual plot since li.pl.dotplot has issues
print("\nCreating summary plot...")

# Count interactions per cell type pair
pair_counts = liana_results.groupby(['source', 'target']).size().reset_index(name='n_interactions')
pair_counts = pair_counts.sort_values('n_interactions', ascending=False).head(20)

fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.barh(
    pair_counts['source'] + ' → ' + pair_counts['target'],
    pair_counts['n_interactions']
)
ax.set_xlabel('Number of interactions')
ax.set_title('Top 20 Cell-Cell Communication Pairs')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/10_liana_summary.png", dpi=150)
print(f"Saved: {FIGURES_DIR}/10_liana_summary.png")
plt.close()

# Save
adata_filtered.write_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_liana.h5ad")
print(f"Saved: {PROCESSED_DIR}/GSE120575_melanoma_liana.h5ad")

print("\nDone!")
