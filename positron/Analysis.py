import pandas as pd
import math
import numpy as np
from compton.Analysis import *
from compton.Dataloader import *
from compton.Functions import *
from positron.Functions import *
from cavendish.utils.Functions import *
from compton.Functions import *
from millikan.functions import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot_counts(df, Normalized = False, overlap=False):
    """Plots the counts as a function of angle.

    Args:
        df (pd.DataFrame): dataframe with index as "angle" and columns "counts" and "error".
    """
    angles = df.index
    counts = df["Counts"]
    errors = df["Error"]
    if Normalized:
        counts = df["normalized_counts"]
        errors = df["normalized_error"]
    plt.errorbar(angles, counts, yerr=errors, marker="o")
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Counts")
    plt.title("Counts vs Angle")
    
    if overlap:
        #plot the overlap
        x = np.linspace(-6, 6, 100)
        overlap = [normalized_overlap(angle) for angle in x]
        plt.plot(x, overlap, label="Overlap", color="red")
        plt.legend()
    
    plt.show()
    
    

def normalize_counts(df):
    """Normalizes the counts in a DataFrame to the maximum count.

    Args:
        df (pd.DataFrame): dataframe with columns "angle" and "counts".

    Returns:
        pd.DataFrame: input dataframe with new collumn of normalized counts, and normalized errors
    """
    max_count = df["Counts"].max()
    print(max_count)
    df["normalized_counts"] = df["Counts"] / max_count
    df["normalized_error"] = df["Error"] / max_count
    return df


def total_counts(dfs):
    #takes in a dictionary of dataframes for each angle, returns a dataframe with the total counts and errors at each angle
    counts = {angle: df["Counts"].sum() for angle, df in dfs.items()} #create a dictionary from the dfs
    counts = pd.DataFrame.from_dict(counts, orient='index', columns=["Counts"]) # create a dataframe from the dictionary
    counts.index.name = "Angle" # set the index name to "Angle"

    counts["Error"] = np.sqrt(counts["Counts"])
    counts = normalize_counts(counts)
    return counts


def run_all_peak_fits(dfs):
    #runs fits for all ROI peaks from Atotal dataframe, returns dataframe with statistics for the peaks
    #args: dfs. this is a DICTIONARY of trimmed dataframes for each  ROI of the form {"ROI": dataframe}
    rows= []
    for roi, df in dfs.items():
        popt, pcov = guassian_fit(df, p0_overide=([175,490,20,] if roi == "High" else None))
        mean, sigma = popt[1], popt[2]
        Unc_mean, Unc_sigma = np.sqrt(np.diag(pcov)[1:3])
        rows.append({"ROI": roi, "Mean": mean, "Sigma": sigma, "Unc Mean": Unc_mean, "Unc Sigma": Unc_sigma})
    peaks = pd.DataFrame(rows)
    peaks.set_index("ROI", inplace=True)
    return peaks