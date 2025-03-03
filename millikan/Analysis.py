"""A Python script containing functions to allow for the SECOND stage of analysis in the millikan oil drop experiment.
The functions here start from a point of having arrays of charges and the uncertainty in said charges, and will include functions to 
plot ideograms, cluster charges based on those ideograms, and plot the clusters in a linear fit, using a weighted average
of the charges in the cluster

***ths script may also include a packaging function to pull together disperate csv's with charge data, and report back as one N*2 array. 
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.signal import find_peaks
import millikan.Dataloader as dl
from millikan.functions import *

def ideogram(q, sigma_q):
    """creates an ideogram for an array of charges and their uncertainties.

    Args:
        q (list): list of droplet charges
        sigma_q (list): corresponding list of droplet charge uncertainties_
    """
    range_q = np.max(q) - np.min(q)
    buffer = 0.2 #fraction of range to use as buffer for plot
    x = np.linspace(np.min(q)-buffer*range_q, np.max(q), 10000)
    y = np.sum([norm.pdf(x, v, u) for v, u in zip(q, sigma_q)], axis=0)
    minima_indices, _ = find_peaks(-y)
    minima_x = x[minima_indices]
    # Plot the summed Gaussians
    plt.plot(x, y, label="Summed Gaussians")
    for min_x in minima_x:
        plt.axvline(x=min_x, color='red', linestyle=':')
    plt.xlabel("charge")
    plt.ylabel("Density")
    plt.title("Ideogram (Summed Gaussian Charge Distributions)")
    plt.legend()
    plt.show()
    
    return minima_x


def bin_charges(q, sigma_q, bin_edges):
    """
    Bins charge values into custom-defined clusters and prepares the data for weighted averaging.

    Args:
        q (array-like): List or array of charge values.
        sigma_q (array-like): Corresponding uncertainties in charge.
        bin_edges (array-like): List of bin boundaries.

    Returns:
        DataFrame: A pandas DataFrame with charge values, their uncertainties, and bin assignments.
    """
    q = np.array(q)
    sigma_q = np.array(sigma_q)

    # Assign each charge value to a bin based on the bin edges
    bin_indices = np.digitize(q, bins=bin_edges, right=False)  

    # Create a DataFrame for easy weighted averaging later
    df = pd.DataFrame({'charge': q, 'sigma_charge': sigma_q, 'bin': bin_indices})

    return df

def weighted_average(df):
    """Computes the weighted average of charge in each bin, weighted by 1/sigma^2."""
    def w_avg(group):
        weights = 1 / group["sigma_charge"]**2
        weighted_mean = np.sum(group["charge"] * weights) / np.sum(weights)
        weighted_sigma = np.sqrt(1 / np.sum(weights))  # Standard error
        return pd.Series({'mean_charge': weighted_mean, 'sigma_mean': weighted_sigma})

    return df.groupby("bin").apply(w_avg).reset_index()