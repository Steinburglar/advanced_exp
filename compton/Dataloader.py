import pandas as pd
import math
import numpy as np
from cavendish.utils.Functions import *
from millikan.functions import *

def load_gamma_count(csv_path):
    """
    Loads gamma count data from CSV file into a nested dictionary.
    
    Args:
        csv_path (str): Path to the CSV file.
        
    Returns:
        dict: Nested dictionary with Channel, Energy, and Counts
    """
    
    # df = pd.read_csv(csv_path, header=6, skiprows=[0,1,2,3,4,5], usecols=[0,1,2])
    
    df = pd.read_csv(csv_path, skiprows=6, names=["Channel", "Energy (keV)", "Counts"])
    
    return df