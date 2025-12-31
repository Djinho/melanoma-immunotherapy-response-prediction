# Predicting Melanoma Immunotherapy Response from Cell-Cell Communication Patterns

[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This project explores whether **cell-cell communication patterns** in the tumour microenvironment can predict response to checkpoint inhibitor immunotherapy in melanoma patients — **before treatment begins**.

While previous studies have characterised cell type differences between responders and non-responders, this analysis focuses on **how immune cells communicate with each other**, hypothesising that the "conversation patterns" between cells may be predictive biomarkers.

### Key Question

> Can pre-treatment cell-cell communication signatures in the tumour microenvironment predict which melanoma patients will respond to anti-PD1 immunotherapy?

---

## Key Findings

### 1. Differential Communication Patterns Identified

Using LIANA (Ligand-Receptor Analysis), we identified communication patterns that differ between responders and non-responders:

| Stronger in Responders | Stronger in Non-Responders |
|------------------------|---------------------------|
| APOE → TREM2 (Macrophages) | IL16 → CCR5 (B cells → T cells) |
| TNF → NOTCH1 (pDC → Monocytes) | IL16 → KCNA3 (B cells → T cells) |
| TGFB1 → ENG (DC2 → Monocytes) | ADAM10 → NOTCH1 (Plasma → B cells) |
| CCL5 → SDC4 (B cells → DC2) | TNF → NOTCH1 (DC2 → B cells) |

**Biological interpretation:** 
- Responders show more activated macrophage and dendritic cell communication
- Non-responders show elevated B cell → T cell immunosuppressive signalling (IL16)

### 2. Regulatory T Cell Communication is Predictive

Machine learning analysis revealed that **Treg ↔ Cytotoxic T cell** communication patterns are among the most important features for predicting response:

- Treg → Cytotoxic T cells (HLA-DRB5 → LAG3)
- Treg → Cytotoxic T cells (CD58 → CD2)
- Treg → Cytotoxic T cells (ADAM10 → NOTCH2)

This suggests Tregs may be actively suppressing anti-tumour T cell responses in non-responders.

### 3. Model Performance

| Model | Accuracy | Features |
|-------|----------|----------|
| Cell type + Gene expression | 78.9% | 34 |
| Communication patterns (LIANA) | 63.2% | 1,281 |

The communication-based model's lower accuracy reflects the challenge of using many features with few samples (19 patients), not necessarily that communication is uninformative.

---

## Dataset

**GSE120575** — Sade-Feldman et al., Cell 2018

> "Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma"

| Metric | Value |
|--------|-------|
| Total cells | 16,291 |
| Pre-treatment cells | 5,927 |
| Patients (pre-treatment) | 19 |
| Responders | 9 |
| Non-responders | 10 |
| Cell types identified | 21 |

---

## Methods

### Pipeline Overview

```
Raw Data (GEO)
    │
    ▼
Quality Control & Filtering
    │
    ▼
Normalisation & HVG Selection
    │
    ▼
Dimensionality Reduction (PCA → UMAP)
    │
    ▼
Clustering (Leiden)
    │
    ▼
Cell Type Annotation (CellTypist)
    │
    ▼
Subset to Pre-treatment Samples
    │
    ▼
Cell-Cell Communication Analysis (LIANA)
    │
    ├──► Group Comparison (Responders vs Non-responders)
    │
    └──► Machine Learning Prediction (Random Forest)
```

### Tools Used

| Tool | Purpose |
|------|---------|
| Scanpy | Single-cell analysis pipeline |
| CellTypist | Automated cell type annotation |
| LIANA | Ligand-receptor communication inference |
| scikit-learn | Machine learning (Random Forest) |

---

## Repository Structure

```
melanoma-immunotherapy-response-prediction/
│
├── README.md
├── data/
│   ├── raw/                    # Original GEO files (not tracked)
│   └── processed/              # Processed h5ad files
│
├── scripts/
│   ├── 01_eda_raw_data.py
│   ├── 02_convert_to_h5ad.py
│   ├── 03_qc_filter.py
│   ├── 04_normalize_hvg.py
│   ├── 05_pca_umap.py
│   ├── 06_clustering.py
│   ├── 07_cell_type_markers.py
│   ├── 07b_celltypist_annotation.py
│   ├── 08_subset_pretreatment.py
│   ├── 09_liana_communication.py
│   ├── 10_compare_responders.py
│   └── 11_ml_prediction.py
│
├── figures/
│   ├── 01_qc_metrics.png
│   ├── 02_highly_variable_genes.png
│   ├── 03_pca_variance.png
│   ├── 04_umap_metadata.png
│   ├── 05_clusters.png
│   ├── 09_celltypist_annotation.png
│   ├── 10_liana_summary.png
│   ├── 11_differential_interactions.png
│   └── 12_ml_communication_results.png
│
├── results/
│   ├── liana_all_interactions.csv
│   ├── liana_responder_vs_nonresponder.csv
│   ├── communication_feature_importance.csv
│   └── ml_predictions.csv
│
└── logs/
    └── *.log
```

---

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/melanoma-immunotherapy-response-prediction.git
cd melanoma-immunotherapy-response-prediction

# Create conda environment
conda create -n melanoma_env python=3.10
conda activate melanoma_env

# Install dependencies
pip install scanpy anndata liana celltypist scikit-learn matplotlib seaborn pandas numpy
pip install igraph leidenalg scikit-misc
```

---

## Usage

Run scripts sequentially:

```bash
cd scripts

# Data processing
python 02_convert_to_h5ad.py
python 03_qc_filter.py
python 04_normalize_hvg.py
python 05_pca_umap.py
python 06_clustering.py

# Cell type annotation
python 07b_celltypist_annotation.py

# Communication analysis
python 08_subset_pretreatment.py
python 09_liana_communication.py
python 10_compare_responders.py

# Machine learning
python 11_ml_prediction.py
```

---

## Limitations

| Limitation | Impact |
|------------|--------|
| **Small sample size (19 patients)** | Insufficient for robust ML conclusions |
| **Single dataset** | No external validation |
| **Communication inferred, not measured** | LIANA infers from gene expression, not direct protein interaction |
| **Class imbalance** | Slightly more non-responders |
| **High dimensionality** | 1,281 features for 19 samples risks overfitting |

---

## Future Directions

### 1. External Validation
Test the identified communication signatures on independent melanoma immunotherapy cohorts (e.g., GSE91061, GSE115978).

### 2. Feature Selection
Reduce features from 1,281 to top 50-100 most differential interactions to improve model stability.

### 3. Multi-cohort Integration
Combine multiple single-cell datasets using batch correction methods (Harmony, scVI) to increase sample size.

### 4. Spatial Validation
Validate key interactions (e.g., Treg → CD8 T cell) using spatial transcriptomics to confirm cells are physically proximal.

### 5. Functional Validation
Experimentally validate IL16 as an immunosuppressive signal in melanoma using in vitro co-culture systems.

---

## Key Conclusions

1. **Cell-cell communication patterns differ** between immunotherapy responders and non-responders at baseline

2. **B cell IL16 signalling** to T cells is elevated in non-responders — a potentially novel immunosuppressive mechanism

3. **Treg ↔ Cytotoxic T cell communication** is among the most predictive features, suggesting active immune suppression in non-responders

4. **Communication-based prediction is feasible** but requires larger cohorts for robust performance

5. **This approach complements** traditional cell type proportion analysis by capturing functional immune interactions

---

## References

1. Sade-Feldman, M. et al. (2018). Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma. *Cell*, 175(4), 998-1013.

2. Dimitrov, D. et al. (2022). Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nature Communications*, 13, 3224.

3. Domínguez Conde, C. et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*, 376(6594).

---

## Author

**Djinho Itshary**

MSc Bioinformatics | Data Scientist

[LinkedIn](https://www.linkedin.com/in/djinho-itshary-671658254/) | [GitHub](https://github.com/Djinho?tab=repositories)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Acknowledgements

- Data from GEO (GSE120575)
- CellTypist for automated cell type annotation
- LIANA developers for the communication inference framework
