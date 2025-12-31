"""
02_convert_to_h5ad.py
Convert raw GSE120575 data to AnnData h5ad format
Memory-efficient version for 16GB RAM systems

Author: Djinho
Date: 2025-12-28
"""

import os
import sys
import logging
import gzip
import gc
from datetime import datetime

import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse

# =============================================================================
# SETUP LOGGING
# =============================================================================
log_dir = "../logs"
os.makedirs(log_dir, exist_ok=True)

log_filename = f"{log_dir}/02_convert_to_h5ad_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# =============================================================================
# CONFIGURATION
# =============================================================================
RAW_DIR = "../data/raw"
PROCESSED_DIR = "../data/processed"
EXPR_FILE = "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
META_FILE = "GSE120575_patient_ID_single_cells.txt.gz"
OUTPUT_FILE = "GSE120575_melanoma.h5ad"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def check_file_exists(filepath, description):
    """Check if file exists and log result."""
    if os.path.exists(filepath):
        size_mb = os.path.getsize(filepath) / (1024 * 1024)
        logger.info(f"✓ Found {description}: {filepath} ({size_mb:.1f} MB)")
        return True
    else:
        logger.error(f"✗ Missing {description}: {filepath}")
        return False


def load_metadata(filepath):
    """Load and clean metadata file."""
    logger.info("Loading metadata...")
    
    try:
        meta = pd.read_csv(filepath, sep='\t', skiprows=19, encoding='latin-1')
        logger.info(f"  Raw metadata shape: {meta.shape}")
        
        cols_to_keep = [
            'title',
            'characteristics: patinet ID (Pre=baseline; Post= on treatment)',
            'characteristics: response',
            'characteristics: therapy'
        ]
        
        meta = meta[cols_to_keep].copy()
        meta.columns = ['cell_id', 'patient_timepoint', 'response', 'therapy']
        
        # Clean up
        meta = meta[meta['response'].isin(['Responder', 'Non-responder'])].copy()
        
        # Extract timepoint and patient
        meta['timepoint'] = meta['patient_timepoint'].apply(
            lambda x: 'Pre' if str(x).startswith('Pre') else 'Post'
        )
        meta['patient_id'] = meta['patient_timepoint'].apply(
            lambda x: str(x).split('_')[-1] if pd.notna(x) else None
        )
        
        meta = meta.set_index('cell_id')
        
        logger.info(f"  Cleaned metadata shape: {meta.shape}")
        logger.info(f"  Response: {meta['response'].value_counts().to_dict()}")
        logger.info(f"  Timepoint: {meta['timepoint'].value_counts().to_dict()}")
        
        return meta
        
    except Exception as e:
        logger.error(f"  Failed to load metadata: {e}")
        raise


def load_expression_memory_efficient(filepath):
    """
    Load expression matrix in a memory-efficient way.
    Reads line by line to avoid loading everything at once.
    """
    logger.info("Loading expression matrix (memory-efficient mode)...")
    logger.info("  This will take several minutes - please be patient...")
    
    try:
        # First pass: get dimensions and cell IDs
        logger.info("  Pass 1: Reading header to get cell IDs...")
        with gzip.open(filepath, 'rt', encoding='latin-1') as f:
            # First line has cell IDs
            header_line = f.readline().strip().split('\t')
            cell_ids = header_line[1:]  # Skip first column (gene name header)
            
            # Second line might be sample IDs - check
            second_line = f.readline().strip().split('\t')
            if second_line[0] == '' or not second_line[1].replace('.', '').replace('-', '').isdigit():
                # Second line is sample info, not expression data
                logger.info("  Detected sample ID row, will skip it...")
                skip_second = True
                sample_ids = second_line[1:]
            else:
                skip_second = False
            
            # Count genes
            n_genes = sum(1 for _ in f)
            if not skip_second:
                n_genes += 1  # Add back the second line we read
                
        n_cells = len(cell_ids)
        logger.info(f"  Found {n_cells} cells, {n_genes} genes")
        
        # Second pass: read expression values
        logger.info("  Pass 2: Reading expression values...")
        
        # Pre-allocate array (float32 to save memory)
        expr_matrix = np.zeros((n_genes, n_cells), dtype=np.float32)
        gene_names = []
        
        with gzip.open(filepath, 'rt', encoding='latin-1') as f:
            # Skip header
            f.readline()
            if skip_second:
                f.readline()
            
            for i, line in enumerate(f):
                if i % 5000 == 0:
                    logger.info(f"    Processing gene {i}/{n_genes}...")
                
                parts = line.strip().split('\t')
                gene_name = parts[0]
                gene_names.append(gene_name)
                
                # Parse expression values
                values = []
                for v in parts[1:n_cells+1]:
                    try:
                        values.append(float(v))
                    except:
                        values.append(0.0)
                
                expr_matrix[i, :len(values)] = values
        
        logger.info(f"  Loaded {len(gene_names)} genes")
        
        # Transpose to cells x genes
        logger.info("  Transposing to cells × genes format...")
        expr_matrix = expr_matrix.T
        
        # Create DataFrame
        logger.info("  Creating DataFrame...")
        expr_df = pd.DataFrame(
            expr_matrix,
            index=cell_ids,
            columns=gene_names
        )
        
        # Free memory
        del expr_matrix
        gc.collect()
        
        logger.info(f"  Final shape: {expr_df.shape}")
        
        return expr_df
        
    except Exception as e:
        logger.error(f"  Failed to load expression matrix: {e}")
        raise


def create_anndata(expr, meta):
    """Create AnnData object from expression and metadata."""
    logger.info("Creating AnnData object...")
    
    try:
        # Find common cells
        common_cells = expr.index.intersection(meta.index)
        logger.info(f"  Cells in expression: {len(expr.index)}")
        logger.info(f"  Cells in metadata: {len(meta.index)}")
        logger.info(f"  Common cells: {len(common_cells)}")
        
        if len(common_cells) == 0:
            logger.error("  No common cells found!")
            logger.error(f"  Expression index examples: {expr.index[:5].tolist()}")
            logger.error(f"  Metadata index examples: {meta.index[:5].tolist()}")
            raise ValueError("No common cells")
        
        # Subset to common cells
        expr_subset = expr.loc[common_cells]
        meta_subset = meta.loc[common_cells]
        
        # Convert to sparse matrix to save memory
        logger.info("  Converting to sparse matrix...")
        X_sparse = sparse.csr_matrix(expr_subset.values)
        
        # Create AnnData
        adata = ad.AnnData(
            X=X_sparse,
            obs=meta_subset,
            var=pd.DataFrame(index=expr_subset.columns)
        )
        
        # Add metadata
        adata.uns['dataset'] = 'GSE120575'
        adata.uns['paper'] = 'Sade-Feldman et al. Cell 2018'
        adata.uns['date_created'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        
        logger.info(f"  AnnData created: {adata.shape[0]} cells × {adata.shape[1]} genes")
        
        # Free memory
        del expr_subset
        gc.collect()
        
        return adata
        
    except Exception as e:
        logger.error(f"  Failed to create AnnData: {e}")
        raise


def save_anndata(adata, filepath):
    """Save AnnData to h5ad file."""
    logger.info(f"Saving to {filepath}...")
    
    try:
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        adata.write_h5ad(filepath)
        
        size_mb = os.path.getsize(filepath) / (1024 * 1024)
        logger.info(f"  ✓ Saved: {filepath} ({size_mb:.1f} MB)")
        
    except Exception as e:
        logger.error(f"  Failed to save: {e}")
        raise


# =============================================================================
# MAIN
# =============================================================================
def main():
    logger.info("=" * 70)
    logger.info("STARTING DATA CONVERSION (Memory-Efficient Mode)")
    logger.info("=" * 70)
    
    # Step 1: Check files
    logger.info("\n[Step 1/5] Checking input files...")
    expr_path = os.path.join(RAW_DIR, EXPR_FILE)
    meta_path = os.path.join(RAW_DIR, META_FILE)
    
    if not (check_file_exists(expr_path, "Expression matrix") and 
            check_file_exists(meta_path, "Metadata file")):
        sys.exit(1)
    
    # Step 2: Load metadata
    logger.info("\n[Step 2/5] Loading metadata...")
    meta = load_metadata(meta_path)
    
    # Step 3: Load expression (memory efficient)
    logger.info("\n[Step 3/5] Loading expression matrix...")
    expr = load_expression_memory_efficient(expr_path)
    
    # Step 4: Create AnnData
    logger.info("\n[Step 4/5] Creating AnnData...")
    adata = create_anndata(expr, meta)
    
    # Free expression dataframe
    del expr
    gc.collect()
    
    # Step 5: Save
    logger.info("\n[Step 5/5] Saving h5ad file...")
    output_path = os.path.join(PROCESSED_DIR, OUTPUT_FILE)
    save_anndata(adata, output_path)
    
    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("CONVERSION COMPLETE!")
    logger.info("=" * 70)
    logger.info(f"Output: {output_path}")
    logger.info(f"Cells: {adata.shape[0]}")
    logger.info(f"Genes: {adata.shape[1]}")
    logger.info(f"Responders: {(adata.obs['response'] == 'Responder').sum()}")
    logger.info(f"Non-responders: {(adata.obs['response'] == 'Non-responder').sum()}")
    logger.info(f"Pre-treatment: {(adata.obs['timepoint'] == 'Pre').sum()}")
    logger.info(f"Post-treatment: {(adata.obs['timepoint'] == 'Post').sum()}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"\nFATAL ERROR: {e}")
        sys.exit(1)
