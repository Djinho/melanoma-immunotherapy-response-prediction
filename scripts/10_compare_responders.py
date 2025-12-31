"""
10_compare_responders.py
Compare cell-cell communication between responders and non-responders
THIS IS THE NOVEL ANALYSIS
"""

import scanpy as sc
import liana as li
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings

warnings.filterwarnings('ignore')

PROCESSED_DIR = "../data/processed"
FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

# Load pre-treatment data
print("Loading pre-treatment data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_pre_treatment.h5ad")
adata.obs['cell_type'] = adata.obs['majority_voting']

# Filter cell types
cell_counts = adata.obs['cell_type'].value_counts()
valid_types = cell_counts[cell_counts >= 30].index.tolist()
adata = adata[adata.obs['cell_type'].isin(valid_types)].copy()

print(f"Cells: {adata.shape[0]}")
print(f"Responders: {(adata.obs['response'] == 'Responder').sum()}")
print(f"Non-responders: {(adata.obs['response'] == 'Non-responder').sum()}")

# Split by response
adata_resp = adata[adata.obs['response'] == 'Responder'].copy()
adata_nonresp = adata[adata.obs['response'] == 'Non-responder'].copy()

# Run LIANA on each group separately
print("\n" + "="*50)
print("Running LIANA on RESPONDERS...")
print("="*50)
li.mt.rank_aggregate(
    adata_resp,
    groupby='cell_type',
    resource_name='consensus',
    use_raw=False,
    verbose=True
)
resp_results = adata_resp.uns['liana_res'].copy()
resp_results['group'] = 'Responder'

print("\n" + "="*50)
print("Running LIANA on NON-RESPONDERS...")
print("="*50)
li.mt.rank_aggregate(
    adata_nonresp,
    groupby='cell_type',
    resource_name='consensus',
    use_raw=False,
    verbose=True
)
nonresp_results = adata_nonresp.uns['liana_res'].copy()
nonresp_results['group'] = 'Non-responder'

# Combine results
print("\nComparing groups...")
all_results = pd.concat([resp_results, nonresp_results], ignore_index=True)

# Create interaction ID for matching
all_results['interaction'] = (
    all_results['source'] + '|' + 
    all_results['target'] + '|' + 
    all_results['ligand_complex'] + '|' + 
    all_results['receptor_complex']
)

# Pivot to compare
pivot = all_results.pivot_table(
    index='interaction',
    columns='group',
    values='magnitude_rank',
    aggfunc='mean'
).dropna()

print(f"Interactions found in both groups: {len(pivot)}")

# Calculate difference (lower rank = stronger interaction)
# Negative diff = stronger in responders
pivot['diff'] = pivot['Non-responder'] - pivot['Responder']
pivot['abs_diff'] = pivot['diff'].abs()

# Top interactions STRONGER in responders (positive diff)
print("\n" + "="*50)
print("TOP 15 INTERACTIONS STRONGER IN RESPONDERS")
print("="*50)
stronger_resp = pivot.sort_values('diff', ascending=False).head(15)
for idx in stronger_resp.index:
    parts = idx.split('|')
    diff = stronger_resp.loc[idx, 'diff']
    print(f"{parts[0]:30} → {parts[1]:30} | {parts[2]:10} → {parts[3]:15} | diff: {diff:.4f}")

# Top interactions STRONGER in non-responders (negative diff)
print("\n" + "="*50)
print("TOP 15 INTERACTIONS STRONGER IN NON-RESPONDERS")
print("="*50)
stronger_nonresp = pivot.sort_values('diff', ascending=True).head(15)
for idx in stronger_nonresp.index:
    parts = idx.split('|')
    diff = stronger_nonresp.loc[idx, 'diff']
    print(f"{parts[0]:30} → {parts[1]:30} | {parts[2]:10} → {parts[3]:15} | diff: {diff:.4f}")

# Save differential results
pivot_save = pivot.reset_index()
pivot_save[['source', 'target', 'ligand', 'receptor']] = pivot_save['interaction'].str.split('|', expand=True)
pivot_save.to_csv(f"{RESULTS_DIR}/liana_responder_vs_nonresponder.csv", index=False)
print(f"\nSaved: {RESULTS_DIR}/liana_responder_vs_nonresponder.csv")

# Plot top differential interactions
print("\nPlotting...")
fig, axes = plt.subplots(1, 2, figsize=(14, 8))

# Stronger in responders
top_resp = pivot.sort_values('diff', ascending=False).head(10)
labels_resp = [f"{x.split('|')[0][:15]}→{x.split('|')[1][:15]}\n({x.split('|')[2]})" for x in top_resp.index]
axes[0].barh(range(len(top_resp)), top_resp['diff'], color='green', alpha=0.7)
axes[0].set_yticks(range(len(top_resp)))
axes[0].set_yticklabels(labels_resp)
axes[0].set_xlabel('Rank difference (+ = stronger in responders)')
axes[0].set_title('Interactions STRONGER in RESPONDERS')
axes[0].invert_yaxis()

# Stronger in non-responders
top_nonresp = pivot.sort_values('diff', ascending=True).head(10)
labels_nonresp = [f"{x.split('|')[0][:15]}→{x.split('|')[1][:15]}\n({x.split('|')[2]})" for x in top_nonresp.index]
axes[1].barh(range(len(top_nonresp)), -top_nonresp['diff'], color='red', alpha=0.7)
axes[1].set_yticks(range(len(top_nonresp)))
axes[1].set_yticklabels(labels_nonresp)
axes[1].set_xlabel('Rank difference (+ = stronger in non-responders)')
axes[1].set_title('Interactions STRONGER in NON-RESPONDERS')
axes[1].invert_yaxis()

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/11_differential_interactions.png", dpi=150, bbox_inches='tight')
print(f"Saved: {FIGURES_DIR}/11_differential_interactions.png")
plt.close()

print("\n" + "="*50)
print("ANALYSIS COMPLETE!")
print("="*50)
