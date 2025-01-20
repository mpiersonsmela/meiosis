import os
from pathlib import Path
import pandas as pd
from dotenv import load_dotenv
import scanpy as sc
import anndata

# Load environment variables
load_dotenv()
DATA_PATH = Path(os.getenv('DATA_PATH'))

def merge_samples():
    # Load reference data
    # Load environment variables
    garcia_path = DATA_PATH / 'garcia_ATAC'
    matrix_files = list(garcia_path.glob('*matrix.mtx.gz'))
     # Load new cell type predictions
    print("Loading cell type predictions...")      
    new_types = pd.read_csv(garcia_path / 'celltype_predictions' / 'all_celltypes.csv', index_col=0)

    # Process each sample
    all_samples = []
    for matrix_file in matrix_files:
        # Get the prefix before _matrix.mtx.gz
        name = matrix_file.stem.replace('matrix.mtx', '')
        # Return None if name contains "and" to skip these samples
        sample_name = None if "and" in name.lower() else name
        if sample_name is None:
            print(f"Skipping {matrix_file.name} as it contains 'and' in the name")
            continue
            
        print(f"Processing sample {sample_name}...")
        
        sample = sc.read_10x_mtx(garcia_path, prefix = sample_name)
        sample.obs["sample"] = sample_name
        new_types_sample = new_types[new_types["sample"] == sample_name]
        sample.obs["celltype"] = new_types_sample["celltype"]
        sample.obs["is_germcell"] = new_types_sample["is_germcell"]
        sample = sample[sample.obs["is_germcell"]]
        all_samples.append(sample)

    all_samples_adata = anndata.concat(all_samples, join = "inner")    

    # Save the result
    print("Saving combined dataset...")
    output_path = garcia_path / 'combined_samples.h5ad'
    all_samples_adata.write_h5ad(output_path)
    print(f"Combined dataset saved to {output_path}")
    
    return all_samples_adata

if __name__ == "__main__":
    merge_samples()
