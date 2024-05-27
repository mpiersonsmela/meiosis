#Code by Laura Breimann with minor changes by Merrick Pierson Smela

# %% Imports
import matplotlib.pyplot as plt
import pandas as pd
import os
import logging
import sys

# %% Configuration
CSV_FOLDER_PATH = sys.argv[1] 
RESULTS_FOLDER_PATH = sys.argv[1]

# %% Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# %% Function to create scatter plots
def create_scatter_plots(csv_file: str, output_folder: str) -> None:
    """Generate and save scatter plots from a CSV file for each fluorescence channel."""
    try:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file)

        # Define colors for the plots
        colors = ['green', 'magenta', 'gray', 'red', 'blue']

        # Plotting
        fig, axs = plt.subplots(1, 5, figsize=(20, 4))  # Adjust the size as necessary
        fig.suptitle(f'Scatter Plots for {os.path.basename(csv_file)}')

        # Loop through each channel to create scatter plots
        for i in range(5):  # Assuming there are 5 channels
            axs[i].scatter(df['area_pixels'], df[f'integrated_intensity_channel_{i}'], color=colors[i])
            axs[i].set_title(f'Channel {i}')
            axs[i].set_xlabel('Area (pixels)')
            axs[i].set_ylabel('Fluorescence Intensity')

        # Tight layout to ensure no overlap of subplots
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        # Save the figure
        fig_filename = f'{os.path.splitext(os.path.basename(csv_file))[0]}_scatter_plots.png'
        fig.savefig(os.path.join(output_folder, fig_filename))
        plt.close(fig)
        logging.info(f"Scatter plots saved to {fig_filename}")

    except Exception as e:
        logging.error(f"Error creating scatter plots for {csv_file}: {e}")

# %% Main execution block
if __name__ == '__main__':
    # Ensure results folder exists
    if not os.path.exists(RESULTS_FOLDER_PATH):
        os.makedirs(RESULTS_FOLDER_PATH)

    # Find CSV files and create plots
    for file in os.listdir(CSV_FOLDER_PATH):
        if file.endswith('.csv'):
            create_scatter_plots(os.path.join(CSV_FOLDER_PATH, file), RESULTS_FOLDER_PATH)

# %%
