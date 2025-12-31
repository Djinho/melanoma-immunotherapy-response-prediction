"""
01_eda_raw_data.py
Exploratory Data Analysis of GSE120575 melanoma dataset
"""

import pandas as pd
import gzip

# =============================================================================
# LOAD METADATA
# =============================================================================
print("=" * 70)
print("LOADING METADATA")
print("=" * 70)

meta_path = "../data/raw/GSE120575_patient_ID_single_cells.txt.gz"

# Read with latin-1 encoding to handle special characters
meta = pd.read_csv(meta_path, sep='\t', skiprows=19, encoding='latin-1')

print(f"\nMetadata shape: {meta.shape}")
print(f"  - {meta.shape[0]} cells")
print(f"  - {meta.shape[1]} columns")

print(f"\nColumn names:")
for i, col in enumerate(meta.columns):
    print(f"  {i}: {col}")

print(f"\nFirst 5 rows:")
print(meta.head())

# =============================================================================
# EXPLORE KEY VARIABLES
# =============================================================================
print("\n" + "=" * 70)
print("KEY VARIABLES FOR OUR ANALYSIS")
print("=" * 70)

# Find response column (look for columns with few unique values)
print("\nSearching for response/treatment info...")
for col in meta.columns:
    n_unique = meta[col].nunique()
    if n_unique < 30:  # Likely categorical
        print(f"\n'{col}' ({n_unique} unique values):")
        print(meta[col].value_counts())

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY FOR PROJECT")
print("=" * 70)

# Try to identify response column
response_col = [c for c in meta.columns if 'response' in c.lower()]
if response_col:
    print(f"\nResponse distribution:")
    print(meta[response_col[0]].value_counts())

# Count unique patients
patient_cols = [c for c in meta.columns if 'patient' in c.lower() or 'patien' in c.lower()]
print(f"\nPatient-related columns: {patient_cols}")

print("\n" + "=" * 70)
print("EDA COMPLETE - Check output above")
print("=" * 70)
