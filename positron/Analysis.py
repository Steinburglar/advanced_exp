import pandas as pd
import math
import numpy as np
from cavendish.utils.Functions import *
from compton.Functions import *
from millikan.functions import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot_counts(counts):
    """Plots the counts as a function of angle.

    Args:
        counts (dict): dictionary of angle (int): counts (int) pairs.
    """
    angles = list(counts.keys())
    counts = list(counts.values())
    plt.errorbar(angles, counts, yerr=np.sqrt(counts), marker = "o")
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Counts")
    plt.title("Counts vs Angle")
    plt.show()