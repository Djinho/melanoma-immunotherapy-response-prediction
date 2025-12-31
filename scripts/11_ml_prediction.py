"""
11_ml_prediction_proper.py
ML model using ACTUAL cell-cell communication scores per patient
(The proper way - no shortcuts)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import liana as li
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
import warnings

warnings.filterwarnings('ignore')

PROCESSED_DIR = "../data/processed"
RESULTS_DIR = "../results"
FIGURES_DIR = "../figures"

# Load pre-treatment data
print("Loading data...")
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE120575_melanoma_pre_treatment.h5ad")
adata.obs['cell_type'] = adata.obs['majority_voting']

# Filter cell types with enough cells globally
cell_counts = adata.obs['cell_type'].value_counts()
valid_types = cell_counts[cell_counts >= 30].index.tolist()
adata = adata[adata.obs['cell_type'].isin(valid_types)].copy()

print(f"Total cells: {adata.shape[0]}")
print(f"Cell types: {len(valid_types)}")

# Get patient list
patients = adata.obs['patient_id'].unique()
print(f"Patients: {len(patients)}")

# Get patient response labels
patient_response = adata.obs.groupby('patient_id')['response'].first().to_dict()

# Define key interactions to track (from our differential analysis)
key_interactions = [
    ('Macrophages', 'Macrophages', 'APOE', 'TREM2'),
    ('pDC', 'Non-classical monocytes', 'TNF', 'NOTCH1'),
    ('DC2', 'Classical monocytes', 'TGFB1', 'ENG'),
    ('Memory B cells', 'DC2', 'CCL5', 'SDC4'),
    ('B cells', 'CD16- NK cells', 'IL16', 'CCR5'),
    ('B cells', 'Tem/Trm cytotoxic T cells', 'IL16', 'CD4'),
    ('B cells', 'Regulatory T cells', 'IL16', 'CCR5'),
    ('Regulatory T cells', 'Tem/Trm cytotoxic T cells', 'TGFB1', 'TGFBR2'),
    ('Macrophages', 'Tem/Trm cytotoxic T cells', 'CD274', 'PDCD1'),  # PD-L1 -> PD-1
]

# Run LIANA per patient
print("\n" + "="*50)
print("Running LIANA PER PATIENT (this will take a while)...")
print("="*50)

patient_scores = {}

for i, patient in enumerate(patients):
    print(f"\nPatient {i+1}/{len(patients)}: {patient}")
    
    # Subset to this patient
    adata_patient = adata[adata.obs['patient_id'] == patient].copy()
    n_cells = adata_patient.shape[0]
    
    # Check cell type distribution for this patient
    ct_counts = adata_patient.obs['cell_type'].value_counts()
    valid_ct = ct_counts[ct_counts >= 5].index.tolist()  # Need at least 5 cells per type
    
    print(f"  Cells: {n_cells}, Valid cell types: {len(valid_ct)}")
    
    if len(valid_ct) < 3:
        print(f"  Skipping - not enough cell types")
        continue
    
    # Filter to valid cell types for this patient
    adata_patient = adata_patient[adata_patient.obs['cell_type'].isin(valid_ct)].copy()
    
    try:
        # Run LIANA
        li.mt.rank_aggregate(
            adata_patient,
            groupby='cell_type',
            resource_name='consensus',
            use_raw=False,
            verbose=False
        )
        
        # Extract results
        results = adata_patient.uns['liana_res']
        
        # Create interaction ID
        results['interaction'] = (
            results['source'] + '|' + 
            results['target'] + '|' + 
            results['ligand_complex'] + '|' + 
            results['receptor_complex']
        )
        
        # Store all interaction scores for this patient
        patient_scores[patient] = results.set_index('interaction')['magnitude_rank'].to_dict()
        print(f"  Interactions found: {len(results)}")
        
    except Exception as e:
        print(f"  Error: {e}")
        continue

print(f"\n\nPatients with LIANA scores: {len(patient_scores)}")

# Build feature matrix from communication scores
print("\nBuilding feature matrix from communication scores...")

# Get all interactions across all patients
all_interactions = set()
for scores in patient_scores.values():
    all_interactions.update(scores.keys())

print(f"Total unique interactions across patients: {len(all_interactions)}")

# Find interactions present in at least 50% of patients
interaction_counts = {}
for interaction in all_interactions:
    count = sum(1 for scores in patient_scores.values() if interaction in scores)
    interaction_counts[interaction] = count

min_patients = len(patient_scores) * 0.5
common_interactions = [k for k, v in interaction_counts.items() if v >= min_patients]
print(f"Interactions in ≥50% of patients: {len(common_interactions)}")

# Build feature matrix
feature_matrix = pd.DataFrame(index=patient_scores.keys(), columns=common_interactions)

for patient, scores in patient_scores.items():
    for interaction in common_interactions:
        if interaction in scores:
            feature_matrix.loc[patient, interaction] = scores[interaction]

# Fill NaN with median (interaction not detected = weak)
feature_matrix = feature_matrix.fillna(feature_matrix.median())

# Convert to numeric
feature_matrix = feature_matrix.astype(float)

print(f"Feature matrix: {feature_matrix.shape[0]} patients × {feature_matrix.shape[1]} interactions")

# Add response labels
feature_matrix['response'] = feature_matrix.index.map(patient_response)
feature_matrix = feature_matrix.dropna(subset=['response'])

# Prepare X and y
X = feature_matrix.drop('response', axis=1)
y = (feature_matrix['response'] == 'Responder').astype(int)

print(f"\nFinal: {X.shape[0]} patients, {X.shape[1]} communication features")
print(f"Responders: {y.sum()}, Non-responders: {len(y) - y.sum()}")

# Scale features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train model
print("\n" + "="*50)
print("Training Random Forest with Leave-One-Out CV...")
print("="*50)

model = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=3)
loo = LeaveOneOut()

y_pred = cross_val_predict(model, X_scaled, y, cv=loo)
y_prob = cross_val_predict(model, X_scaled, y, cv=loo, method='predict_proba')[:, 1]

# Evaluate
accuracy = accuracy_score(y, y_pred)
print(f"\nAccuracy: {accuracy:.1%}")

print("\nClassification Report:")
print(classification_report(y, y_pred, target_names=['Non-responder', 'Responder']))

print("Confusion Matrix:")
cm = confusion_matrix(y, y_pred)
print(cm)

# Feature importance
print("\n" + "="*50)
print("TOP 15 MOST IMPORTANT COMMUNICATION PATTERNS")
print("="*50)

model.fit(X_scaled, y)
importance = pd.DataFrame({
    'interaction': X.columns,
    'importance': model.feature_importances_
}).sort_values('importance', ascending=False)

for i, row in importance.head(15).iterrows():
    parts = row['interaction'].split('|')
    print(f"{parts[0]:25} → {parts[1]:25} ({parts[2]}→{parts[3]}): {row['importance']:.4f}")

# Save
importance.to_csv(f"{RESULTS_DIR}/communication_feature_importance.csv", index=False)

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Confusion matrix
axes[0].imshow(cm, cmap='Blues')
axes[0].set_xticks([0, 1])
axes[0].set_yticks([0, 1])
axes[0].set_xticklabels(['Non-resp', 'Resp'])
axes[0].set_yticklabels(['Non-resp', 'Resp'])
axes[0].set_xlabel('Predicted')
axes[0].set_ylabel('Actual')
axes[0].set_title(f'Confusion Matrix\nAccuracy: {accuracy:.1%}')
for i in range(2):
    for j in range(2):
        axes[0].text(j, i, cm[i, j], ha='center', va='center', fontsize=20)

# ROC curve
fpr, tpr, _ = roc_curve(y, y_prob)
roc_auc = auc(fpr, tpr)
axes[1].plot(fpr, tpr, 'b-', lw=2, label=f'AUC = {roc_auc:.2f}')
axes[1].plot([0, 1], [0, 1], 'k--', lw=1)
axes[1].set_xlabel('False Positive Rate')
axes[1].set_ylabel('True Positive Rate')
axes[1].set_title('ROC Curve')
axes[1].legend()

# Top interactions
top_imp = importance.head(10)
labels = [f"{x.split('|')[0][:12]}→{x.split('|')[1][:12]}" for x in top_imp['interaction']]
axes[2].barh(range(len(top_imp)), top_imp['importance'])
axes[2].set_yticks(range(len(top_imp)))
axes[2].set_yticklabels(labels)
axes[2].set_xlabel('Importance')
axes[2].set_title('Top Communication Features')
axes[2].invert_yaxis()

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/12_ml_communication_results.png", dpi=150, bbox_inches='tight')
print(f"\nSaved: {FIGURES_DIR}/12_ml_communication_results.png")
plt.close()

print("\n" + "="*50)
print("PROPER ML ANALYSIS COMPLETE!")
print("="*50)
print("\nThis model uses ACTUAL cell-cell communication scores,")
print("not just cell type proportions. This is the novel contribution.")
