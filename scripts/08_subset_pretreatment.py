"""
08_subset_pretreatment.py
Subset to pre-treatment samples only (baseline prediction)
"""

import scanpy as sc
import pandas as pd

PROCESSED_DIR = "../data/processed"

# Load CellTypist annotated data
print("Loading CellTypist annotated data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_celltypist.h5ad")
print(f"Total cells: {adata.shape[0]}")

# Check timepoint distribution
print("\nTimepoint distribution:")
print(adata.obs['timepoint'].value_counts())

# Subset to PRE-treatment only
print("\nSubsetting to pre-treatment only...")
adata_pre = adata[adata.obs['timepoint'] == 'Pre'].copy()
print(f"Pre-treatment cells: {adata_pre.shape[0]}")

# Check response distribution in pre-treatment
print("\nResponse distribution (pre-treatment):")
print(adata_pre.obs['response'].value_counts())

# Check cell type distribution
print("\nCell types in pre-treatment:")
print(adata_pre.obs['majority_voting'].value_counts())

# Check cells per patient
print("\nCells per patient:")
patient_counts = adata_pre.obs.groupby(['patient_id', 'response']).size().reset_index(name='n_cells')
print(patient_counts.to_string())

# How many patients per response group?
print("\nPatients per response group:")
patients_per_response = adata_pre.obs.groupby('response')['patient_id'].nunique()
print(patients_per_response)

# Save
output = f"{PROCESSED_DIR}/GSE120575_melanoma_pre_treatment.h5ad"
adata_pre.write_h5ad(output)
print(f"\nSaved: {output}")
