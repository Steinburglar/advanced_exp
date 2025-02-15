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


def load_convert(path, sig_x = 0.0005, num = 4 ):
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
    df = pd.read_csv(path, header=[0,1])
    dfs = clean_split(df, num, sig_x)
    convert_rad(dfs)
    return dfs

def clean_split(df, num, sig_x):
    """takes in a dataframe, cleans it to keep only the Time (sec) and Position (cm) columns, and splits into 4 dataframes, one for each measurement. 
    Args:
        df (_dataframe_): pandas dataframe to be cleaned and split
        num (_int_): number of seperate measurements to split the data into. assumed to be 4
        sig_x (_float_): uncertainty in position measurement for each timepoint, assumed to be uniform across all observations
    returns:
        List of dataframes, one for each measurement. 
    """
    df = df.loc[:, pd.IndexSlice[:, ["Time (sec)", "Position (cm)"]]]
    df =df.dropna()
    dfs = [df.loc[:, pd.IndexSlice[f"Measurement {i}", :]].copy() for i in range(1, num+1)]
    del df #for memory
    for i, _df in enumerate(dfs):
        _df.loc[:, (f"Measurement {i+1}", "Uncertainty (m)")] = sig_x
        _df.rename(columns={"Position (cm)": "Position (m)"}, level = 1, inplace = True)
        _df.loc[:, _df.columns.get_level_values(1) == 'Position (m)'] /= 100    
    return dfs


def convert_rad(dfs):
    """takes in a list of dataframes, and converts the collums of each from position to radians, based off of our measurements

    Args:
        dfs (_list_): list if pandas dataframes, one for each measurement
    """
    for i, df in enumerate(dfs):
        measure = df.columns.get_level_values(0)[0]
        df[[(measure, "Radians"), (measure, "Uncertainty (rad)")]] = \
            df.apply(lambda row: pd.Series(
                con_rad((row[(measure, "Position (m)")],row[(measure, "Uncertainty (m)")])),            \
                    index=[(measure, "Radians"), (measure, "Uncertainty (rad)")]), axis = 1,)
        df.drop([(measure, "Position (m)"), (measure, "Uncertainty (m)")], axis = 1, inplace=True)
    return
