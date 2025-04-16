"""
Lucas Steinberger
Created on 4/13/2025
Cotains dataloading functions for positron annihillations
"""
import pandas as pd
import math
import numpy as np
from cavendish.utils.Functions import *
from millikan.functions import *
from compton.Functions import *
from compton.Analysis import *
from compton.Dataloader import *

def load_positron_angles(rel_path):
    """
    Loads select positron annihilation data from CSV files into a dictionary of dataframes.
    
    Returns:
        list: dictionary of dataframes containing the loaded data.
    """
    neg_trials = [-1, -2, -3, -4, -5, -6]
    trials = [0, 1, 2, 3, 4, 5, 6]
    total_trials = sorted(neg_trials + trials)
    data = {trial: load_gamma_count(rel_path + "A" + str(trial) + ".csv") for trial in total_trials}
    return data
    