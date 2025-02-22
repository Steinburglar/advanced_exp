"""
Lucas Steinberger
1/31/2025
Python script containing  functions to for loading and converting data from Casey, Lucas, and Cooke's Cavendish lab in PHYS 064.
Goal is to provide utility to load csv and return 4 arrays, each containing the data points with corresponding uncertainties for a measurement,
with position converted to radians of beam rotation.
"""
import pandas as pd
import math
import numpy as np
import os
from cavendish.utils.Functions import *


def bbl_loader(directory):
    """
    Loads all CSV files in a given directory, assuming they contain drop velocity data 
    structured with a hierarchical column format (BBL as the first header and drop numbers as the second header).

    Args:
        directory (str): Path to the folder containing CSV files.

    Returns:
        dict: A dictionary where keys are filenames and values are the corresponding DataFrames.
    """
    dataframes = {}

    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            file_path = os.path.join(directory, filename)
            dataframes[filename] = load_convert(file_path)
    
    return dataframes

def load_convert(path):
    """
    Loads a CSV file containing drop velocity measurements across multiple trials (BBL).

    The CSV has a hierarchical column structure:
        - First header row: 'BBL' marking the trial number.
        - Second header row: 'drop # velo' indicating the drop number for velocity measurements.

    This function loads the data into a pandas DataFrame with a MultiIndex for columns, 
    allowing easy access to velocities for each drop within each trial.

    Args:
        path (str): Path to the .csv file.

    Returns:
        pd.DataFrame: A DataFrame with a MultiIndex for columns (BBL, drop # velo).
    """
    df = pd.read_csv(path, header=[0, 1])  # Read the CSV with two header rows
    df.columns = pd.MultiIndex.from_tuples(df.columns)  # Ensure columns are a MultiIndex
    return df