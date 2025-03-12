"""
Lucas Steinberger
3/8/2025
Python script containing  functions to for loading and converting data from Casey, Lucas, and Cooke's Millikan lab in PHYS 064.
There may actually be no data that needs to be loaded here, but I'm not sure yet.

"""
import pandas as pd
import numpy as np
def create_dataframe(csv_path):
    """creates a dataframe based on data that is in a CSV that lucas assembled of the points where the turns.
    returns: a dataframe with 4  columns: 
        img: for the name of the image from which that data point was pulled,
        n: for the integer number of of the "turn"
        V, for the voltage supplied at the "turn",
        sigma_V: for the uncertainty in the voltage.
    """
    df = pd.read_csv(csv_path, header=0)
    df["big_sigma"] = np.ones(len(df))
    scale = df['V_max'] / 50
    df['V'] = df['x'] * scale
    df['sigma_V'] = df['sigma_x'] * scale
    df['big_sigma_v'] = df['big_sigma'] * scale *2
    return df