import os
from pathlib import Path
import scanpy as sc
import anndata
from dotenv import load_dotenv

load_dotenv()
DATA_PATH = Path(os.getenv('DATA_PATH'))

def reduce_dimensions(file_type):
    adata = sc.read_h5ad(DATA_PATH / 'atlas' / f'processed_files/E-MTAB-10551/human_{file_type}.h5ad')
    adata.X = adata.raw.X
    adata.raw = None

    adata.obs = adata.obs[["celltype"]]
    adata.var = adata.var[["gene_ids-0"]]
    sc.write(DATA_PATH / 'atlas' / f'processed_files/E-MTAB-10551/human_{file_type}_reduced.h5ad', adata)
    return

def create_combined_reference():
    """
    Creates a combined reference dataset from germcell and somatic cell data.
    Returns the combined AnnData object and saves it to disk.
    """
    load_dotenv()

    # Load reference datasets
    ref_germcell = sc.read_h5ad(DATA_PATH / 'atlas' / 'processed_files/E-MTAB-10551/human_germcells_reduced.h5ad')
    ref_somatic = sc.read_h5ad(DATA_PATH / 'atlas' / 'processed_files/E-MTAB-10551/human_somatic_reduced.h5ad')

    # Concatenate the datasets
    ref_combined = anndata.concat([ref_somatic, ref_germcell], join="inner")

    # Normalize and process the combined data
    sc.pp.normalize_total(ref_combined)
    sc.pp.log1p(ref_combined)
    sc.pp.pca(ref_combined)
    sc.pp.neighbors(ref_combined)
    sc.tl.umap(ref_combined)

    # Save the combined reference
    output_path = DATA_PATH / 'atlas' / 'processed_files/E-MTAB-10551/human_combined.h5ad'
    sc.write(output_path, ref_combined)

    return ref_combined

if __name__ == "__main__":
    reduce_dimensions("somatic")
    print("Somatic reduced")
    reduce_dimensions("germcells")
    print("Germcells reduced")
    create_combined_reference() 
    print("Combined reference created")