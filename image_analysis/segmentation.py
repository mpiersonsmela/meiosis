#Code by Laura Breimann with minor changes by Merrick Pierson Smela

# %% Imports
import sys
import os
import logging
from typing import Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cellpose import models
import czifile
from skimage import measure


# %% Paths and model
FOLDER_PATH = sys.argv[1]
dapi_channel = int(sys.argv[2]) # Which channel contains DAPI (channel 0 = first channel)
OUTPUT_FOLDER = sys.argv[3] #'/Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Images/Segmentation_analysis/'
MODEL_TYPE = 'nuclei'  


# %% Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# %% Helper functions
def load_image(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load CZI file and return the DAPI channel data."""
    try:
        img = czifile.imread(file_path)
        #dapi_data = img[0, 0, 2, 0, :, :].squeeze()  # Assumes DAPI is the third channel
        dapi_data = img[0, 0, dapi_channel, 0, :, :].squeeze()  
        return img, dapi_data
    except Exception as e:
        logging.error(f"Failed to read {file_path}: {e}")
        raise

def plot_masks_on_image(image: np.ndarray, masks: np.ndarray, filename: str) -> None:
    """Generate and save an image showing the segmentation masks overlaid on the DAPI channel."""
    plt.figure(figsize=(10, 10))
    plt.imshow(image, cmap='gray')
    regions = measure.regionprops(masks)
    for region in regions:
        plt.text(region.bbox[1], region.bbox[0], region.label, color='white', fontsize=12)
        plt.imshow(np.ma.masked_where(masks == 0, masks), alpha=0.5, cmap='jet')
    plt.title(f'Segmentation for {filename}')
    plt.axis('off')
    plt.savefig(os.path.join(OUTPUT_FOLDER, f'{filename}_segmentation.png'))
    plt.close()

def process_files(folder_path: str, output_folder: str):
    """Process each CZI file in the specified folder."""
    for filename in os.listdir(folder_path):
        if filename.endswith('.czi'):
            if os.path.isfile(os.path.join(output_folder, f'{filename}_fluorescence.csv')):
                print("Already done: "+filename)
            else:
                file_path = os.path.join(folder_path, filename)
                img, dapi_data = load_image(file_path)
                masks, flows, styles, diams = MODEL.eval([dapi_data], diameter=150, channels=[0,0]) #set default nuclei diameter to 150 pixels
                plot_masks_on_image(dapi_data, masks[0], filename)
                df = extract_fluorescence_data(img, masks[0])
                df.to_csv(os.path.join(output_folder, f'{filename}_fluorescence.csv'), index=False)

def extract_fluorescence_data(img: np.ndarray, masks: np.ndarray) -> pd.DataFrame:
    """Extract and store fluorescence intensity and area data from all channels."""
    df_list = []
    for c in range(img.shape[2]):  # Consider each channel
        channel_data = img[0, 0, c, 0, :, :].squeeze()
        regions = measure.regionprops(masks, channel_data)
        data = {
            'label': [region.label for region in regions],
            'area_pixels': [region.area for region in regions],
            f'integrated_intensity_channel_{c}': [region.intensity_image.sum() for region in regions]
        }
        df_channel = pd.DataFrame(data)
        if df_list:
            df = df_list[0].merge(df_channel, on=['label', 'area_pixels'], how='outer')
            df_list[0] = df  # Update the first dataframe with merged results
        else:
            df_list.append(df_channel)

    return df_list[0] if df_list else pd.DataFrame()



# %% 
if __name__ == '__main__':
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
    MODEL = models.Cellpose(gpu=True, model_type=MODEL_TYPE)
    process_files(FOLDER_PATH, OUTPUT_FOLDER)

# %%
