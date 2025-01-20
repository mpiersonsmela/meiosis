import pandas as pd
import os
from pathlib import Path
import logging
from tqdm import tqdm
from dotenv import load_dotenv
import ftplib
import socket
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def download_ftp_file(url, output_path):
    """
    Download a file from an FTP URL to the specified output path with a progress bar
    """
    try:
        # Parse FTP URL
        # Remove 'ftp://' and split into server and path
        url = url.replace('ftp://', '')
        server = url.split('/', 1)[0]
        path = '/' + url.split('/', 1)[1]
        
        # Connect to FTP server
        ftp = ftplib.FTP(server, timeout=60)
        ftp.login()  # anonymous login
        
        # Get file size
        ftp.voidcmd('TYPE I')  # binary mode
        size = ftp.size(path)
        
        # Download with progress bar
        with open(output_path, 'wb') as f:
            with tqdm(total=size, unit='B', unit_scale=True, 
                     desc=os.path.basename(output_path)) as pbar:
                
                def callback(data):
                    f.write(data)
                    pbar.update(len(data))
                
                ftp.retrbinary(f'RETR {path}', callback)
        
        ftp.quit()
        return True
        
    except (ftplib.error_perm, socket.timeout, socket.gaierror, Exception) as e:
        logger.error(f"Error downloading {url}: {str(e)}")
        # Clean up partial download if it exists
        if output_path.exists():
            output_path.unlink()
        return False

def get_sample_pairs(sdrf_path):
    """
    Get pairs of RNA and ATAC samples and their FASTQ URLs starting from RNA samples
    """
    # Read SDRF file
    sdrf = pd.read_csv(sdrf_path, sep='\t')
    
    # Create dictionary of RNA-ATAC pairs by looking at RNA samples
    pairs_dict = {}
    rna_samples = sdrf.loc[
        sdrf["Characteristics[protocol]"] == "rna",
        ["Source Name", "Characteristics[paired library]"]
    ].drop_duplicates()
    
    pairs_dict = dict(zip(
        rna_samples["Source Name"],
        rna_samples["Characteristics[paired library]"]
    ))
    
    # Get FASTQ URLs for each sample (both RNA and ATAC)
    sample_urls = {}
    all_samples = set(pairs_dict.keys()) | set(pairs_dict.values())
    for sample in all_samples:
        urls = sdrf.loc[
            sdrf['Source Name'] == sample,
            'Comment[FASTQ_URI]'
        ].dropna().unique()
        sample_urls[sample] = urls
    
    return pairs_dict, sample_urls

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Download RNA and ATAC FASTQ files for a sample pair')
    parser.add_argument('rna_sample', help='RNA sample name (e.g., FCA_GND10287600)')
    args = parser.parse_args()
    
    # Set paths
    load_dotenv()
    data_dir = Path(os.getenv('DATA_PATH')) / 'garcia_ATAC'
    sdrf_path = data_dir / "E-MTAB-11708.sdrf.txt"
    
    # Get sample pairs and their URLs
    pairs_dict, sample_urls = get_sample_pairs(sdrf_path)
    
    # Check if RNA sample exists
    if args.rna_sample not in pairs_dict:
        logger.error(f"RNA sample {args.rna_sample} not found")
        return
    
    # Get ATAC pair
    atac_sample = pairs_dict[args.rna_sample]
    
    # Create directory for the pair and subfolders
    pair_dir = data_dir / f"{args.rna_sample}"
    rna_dir = pair_dir / "rna"
    atac_dir = pair_dir / "atac"
    rna_dir.mkdir(parents=True, exist_ok=True)
    atac_dir.mkdir(parents=True, exist_ok=True)
    
    # Download RNA files
    for url in sample_urls[args.rna_sample]:
        filename = url.split('/')[-1]
        # Replace original sample name with RNA sample name
        filename = filename.replace(args.rna_sample, args.rna_sample)
        output_path = rna_dir / filename
        
        # Skip if file already exists
        if output_path.exists():
            logger.info(f"Skipping {filename} - already exists")
            continue
        
        logger.info(f"Downloading RNA file {filename} for {args.rna_sample}")
        success = download_ftp_file(url, output_path)
        
        if success:
            logger.info(f"Successfully downloaded {filename}")
        else:
            logger.error(f"Failed to download {filename}")
    
    # Download ATAC files
    for url in sample_urls[atac_sample]:
        filename = url.split('/')[-1]
        # Replace ATAC sample name with RNA sample name
        filename = filename.replace(atac_sample, args.rna_sample)
        output_path = atac_dir / filename
        
        # Skip if file already exists
        if output_path.exists():
            logger.info(f"Skipping {filename} - already exists")
            continue
        
        logger.info(f"Downloading ATAC file {filename} for {atac_sample}")
        success = download_ftp_file(url, output_path)
        
        if success:
            logger.info(f"Successfully downloaded {filename}")
        else:
            logger.error(f"Failed to download {filename}")

if __name__ == "__main__":
    main()
