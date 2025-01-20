import os
from pathlib import Path
import pandas as pd
import scanpy as sc
from dotenv import load_dotenv

def load_reference_data(data_path):
    """Load and prepare the reference data"""
    ref_combined = sc.read_h5ad(data_path / 'atlas' / 'processed_files/E-MTAB-10551/human_combined.h5ad')
    return ref_combined

def process_sample(garcia_path, sample_name, ref_data):
    """Process a single sample and return predicted cell types"""
    # Load the sample data using the directory path
    sample = sc.read_10x_mtx(garcia_path, prefix = sample_name)
    
    # Normalize and process the sample
    sc.pp.normalize_total(sample)
    sc.pp.log1p(sample)
    sc.pp.neighbors(sample)
    sc.tl.umap(sample)
    
    # Get common genes between reference and sample
    common_genes = list(set(ref_data.var_names).intersection(set(sample.var_names)))
    
    # Subset both datasets to common genes
    sample = sample[:, common_genes]
    ref_subset = ref_data[:, common_genes]
    
    # Perform cell type prediction
    sc.tl.ingest(sample, ref_subset, obs="celltype", inplace = True)
    
    return pd.Series(sample.obs["celltype"], index=sample.obs_names)

def is_germcell(celltype):
    """Determine if a cell type is a germ cell"""
    germcell_types = {
        'PGC', 'oogonia_STRA8', 'oogonia_meiotic', 'pre_oocyte',
        'oocyte', 'pre_spermatogonia', 'GC_mitotic', 'GC'
    }
    return celltype in germcell_types

def get_sample_name(matrix_file):
    """Extract sample name from matrix file path"""
    # Get the prefix before _matrix.mtx.gz
    name = matrix_file.stem.replace('matrix.mtx', '')
    # Return None if name contains "and" to skip these samples
    return None if "and" in name.lower() else name

def main():
    # Load environment variables
    load_dotenv()
    DATA_PATH = Path(os.getenv('DATA_PATH'))
    
    # Load reference data
    print("Loading reference data...")
    ref_data = load_reference_data(DATA_PATH)
    
    # Get matrix.mtx files directly in garcia_ATAC folder (no subfolders)
    garcia_path = DATA_PATH / 'garcia_ATAC'
    matrix_files = list(garcia_path.glob('*matrix.mtx.gz'))
    
    if not matrix_files:
        print("No matrix.mtx files found in the garcia_ATAC directory!")
        return
    
    # Initialize empty DataFrame to store all results
    all_predictions = pd.DataFrame()
    
    # Process each sample
    for matrix_file in matrix_files:
        sample_name = get_sample_name(matrix_file)
        if sample_name is None:
            print(f"Skipping {matrix_file.name} as it contains 'and' in the name")
            continue
            
        print(f"Processing sample {sample_name}...")
        
        try:
            predictions = process_sample(garcia_path, sample_name, ref_data)
            
            # Add sample name as a column
            sample_predictions = pd.DataFrame({
                'sample': sample_name,
                'celltype': predictions,
                'is_germcell': predictions.map(is_germcell)
            })
            
            all_predictions = pd.concat([all_predictions, sample_predictions])
            
        except Exception as e:
            print(f"Error processing {sample_name}: {str(e)}")
    
    # Filter for germ cells only
    germcell_predictions = all_predictions[all_predictions['is_germcell']]
    
    # Save results
    output_path = garcia_path / 'celltype_predictions'
    output_path.mkdir(exist_ok=True)
    
    # Save all predictions
    all_predictions.to_csv(output_path / 'all_celltypes.csv')
    
    # Save germ cell predictions
    germcell_predictions.to_csv(output_path / 'germcell_types.csv')
    
    print("Processing complete. Results saved to celltype_predictions folder.")

if __name__ == "__main__":
    main() 