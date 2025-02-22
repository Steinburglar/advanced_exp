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
from cavendish.utils.Functions import *

def load_velocity_data(csv_path):
    """
    Loads velocity data from a CSV file into a nested dictionary.
    
    Args:
        csv_path (str): Path to the CSV file.
        
    Returns:
        dict: Nested dictionary with BBL as keys, Drop numbers as subkeys, and lists of velocities as values.
    """
    df = pd.read_csv(csv_path, header=[0, 1])  # Read with two header rows
    
    velocity_data = {}
    
    # Extract unique BBL labels from the first header row
    unique_bbls = df.columns.get_level_values(0).unique()
    
    for bbl in unique_bbls:
        velocity_data[bbl] = {}
        
        # Filter columns belonging to this BBL
        bbl_columns = df[bbl]
        
        for drop in bbl_columns.columns:  # Second-level headers
            velocity_data[bbl][drop] = bbl_columns[drop].dropna().tolist()  # Store non-null values

    return velocity_data

def load_convert(path, num = 11):
    """ Main function for dataloader. cleans unnecessary columns, splits into separate data frames, adds uncertainty, coverts to radians of beam rotation. 
    Note that that means the values returned are *half* of the angle produced by the laser, since the law of reflection
    will double the angle between the laser rays.

    Args:
        path (_str_): Path to the location of the .csv file containing measurements of time and position (x) in cm.
    Kwargs:
        sig_x (_float_): uncertainty in position measurement for each timepoint, assumed to be uniform across all observations.
        num (_int_): number of seperate measurements to split the data into. assumed to be 4
    Returns:
        _list_: list of data arrays:  [data_array1, data_array2, data_array3 ....]
            _data_array_: numpy array containing 3 collumns:
                Time    | angle  | uncertainty in angle
    """
    df = pd.read_csv(path, header=[0,1], index_col="BBL"+str())
    
    return df



