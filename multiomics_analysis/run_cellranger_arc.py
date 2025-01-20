#!/usr/bin/env python3

import os
import argparse
import subprocess
import json
from pathlib import Path
from dotenv import load_dotenv
import pandas as pd

def create_library_file(sample_name, fastq_dir):
    """
    Create a library CSV file for Cell Ranger ARC
    """
    rna_dir = fastq_dir / "rna"
    atac_dir = fastq_dir / "atac"
    
    library_content = [
        "fastqs,sample,library_type",
        f"{rna_dir},{sample_name},Gene Expression",
        f"{atac_dir},{sample_name},Chromatin Accessibility"
    ]
    
    library_path = f"{sample_name}_library.csv"
    with open(library_path, "w") as f:
        f.write("\n".join(library_content))
    
    # Verify directories exist
    if not rna_dir.exists():
        raise ValueError(f"RNA directory does not exist: {rna_dir}")
    if not atac_dir.exists():
        raise ValueError(f"ATAC directory does not exist: {atac_dir}")
    
    return library_path

def run_cellranger_arc(sample_name, reference_path):
    """
    Run Cell Ranger ARC command for the given sample
    """
    # Get data directory from environment variable
    data_dir = Path(os.getenv('DATA_PATH')) / 'garcia_ATAC'
    
    # Construct fastq directory path
    fastq_dir = data_dir / f"{sample_name}"
    
    if not fastq_dir.exists():
        raise ValueError(f"FASTQ directory does not exist: {fastq_dir}")
    
    # Create library file
    library_path = create_library_file(sample_name, fastq_dir)
    
    # Construct Cell Ranger ARC command
    cmd = [
        "cellranger-arc", "count",
        f"--id={sample_name}",
        f"--reference={reference_path}",
        f"--libraries={library_path}",
        "--disable-ui"
    ]
        
    # Run the command
    try:
        subprocess.run(cmd, check=True)
        print(f"Cell Ranger ARC analysis completed successfully for {sample_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Cell Ranger ARC: {e}")
        raise
    finally:
        # Clean up library file
        os.remove(library_path)

def main():
    parser = argparse.ArgumentParser(description="Run Cell Ranger ARC analysis for multiome data")
    parser.add_argument("sample_name", help="Name of the sample to process")
    parser.add_argument("--reference", 
                       help="Path to Cell Ranger ARC reference",
                       default=str(Path.home() / "refdata-cellranger-arc-GRCh38-2020-A-2.0.0"))
    
    args = parser.parse_args()
    
    # Load environment variables
    load_dotenv()
    
    # Verify reference exists
    reference_path = Path(args.reference)
    if not reference_path.exists():
        raise ValueError(f"Reference directory does not exist: {reference_path}\n"
                        f"Please download it from: https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz")
    
    # Run Cell Ranger ARC
    run_cellranger_arc(args.sample_name, reference_path)

if __name__ == "__main__":
    main()
