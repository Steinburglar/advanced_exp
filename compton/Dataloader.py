import pandas as pd
import math
import numpy as np
from compton.Functions import *
from cavendish.utils.Functions import *
from millikan.functions import *

def load_gamma_count(csv_path, rescale=None):
    """
    Loads gamma count data from CSV file into a dataframe
    
    Args:
        csv_path (str): Path to the CSV file.
        
    Returns:
        dataframe: dataframe with Channel, Energy, and Counts
    """
    
    # df = pd.read_csv(csv_path, header=6, skiprows=[0,1,2,3,4,5], usecols=[0,1,2])
    
    df = pd.read_csv(csv_path, skiprows=7, names=["Channel", "Energy (keV)", "Counts"])
    if rescale is not None:
        df["Energy (keV)"] = linear(df["Energy (keV)"], rescale[0], rescale[1])
    
    return df

def expand_df(df, E, sig_E):
    """expands a dataframe to calculate the values of interest for linear regression, namely E/E', corresponding uncertainty in E/E',
    and 1/cos(theta).

    Args:
        df (_type_): Dataframe containing following collumns: Angle, Mean, Sigma, Unc Mean, Unc Sigma
        E (float): Value of peak energy taken by fitting a guassian to the A0_S0 baseline run 
        sig_E (float): Uncertainty in E, as reported by curve_fit
    """
    df["E/E'"] = E/df["Mean"]
    df["Unc E/E'"] = np.sqrt((df["Unc Mean"]/E)**2 + (df["Mean"]*sig_E/(E**2))**2)
    df["1-cos(theta)"] = 1-np.cos(np.radians(df["Angle"]))
    return df