import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_ecg_data_from_folder(folder_path):
    # Get a list of all CSV files in the specified folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    # Regular expression to extract reset time and height from the CSV filenames
    regex_pattern = r'electrode_data_reset_(\d+)_height_(\d+).csv'

    # Create a subplot for all electrodes on a single plot
    fig, ax = plt.subplots(figsize=(12, 6))

    for csv_file in csv_files:
        # Extract reset time and height from the CSV filename using regex
        match = re.match(regex_pattern, csv_file)
        if match:
            # Load ECG data from CSV
            ecg_data = pd.read_csv(os.path.join(folder_path, csv_file), index_col=0)

            # Resetting the index to use row numbers as x-axis values
            ecg_data_reset = ecg_data.reset_index()

            # Plot each electrode's ECG data
            for electrode_label in ecg_data_reset.columns[1:]:
                ax.plot(ecg_data_reset.index, ecg_data_reset[electrode_label], label=f'{electrode_label} ({match.group(1)}, {match.group(2)})')



    plt.xlabel("Time Step")
    plt.ylabel("Integrated Signal")
    plt.title("ECG Signals Over Time - All Electrodes")
    plt.legend()
    plt.show() 

# Example usage with a specific folder path
ecg_csv_folder = '/Users/alexander/Desktop/MEng_Indv_Project/Simulations'
plot_ecg_data_from_folder(ecg_csv_folder)
