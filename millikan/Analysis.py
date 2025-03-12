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
from scipy.optimize import curve_fit
import millikan.Dataloader as dl
from millikan.functions import *

def ideogram(q, sigma_q):
    """creates an ideogram for an array of charges and their uncertainties.

    Args:
        q (np.array): list of droplet charges
        sigma_q (np.array): corresponding list of droplet charge uncertainties_
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


def prune_data(charges, sigma_charges, threshold):
    pruned_charges = []
    pruned_sigma_charges = []

    for charge, sigma_charge in zip(charges, sigma_charges):
        fractional_error = sigma_charge / charge
        if fractional_error <= threshold:
            pruned_charges.append(charge)
            pruned_sigma_charges.append(sigma_charge)

    return np.array(pruned_charges), np.array(pruned_sigma_charges)

def plot_fractional_error(q, sigma_q, threshold=0.06):
    """plots fractional error in charge as a function of charge value on a log-linear scale.
        it also plots a horizontal line at 1% error for reference, and to help separate the data into two categories.
    Args:
        q (_type_): _description_
        sigma_q (_type_): _description_
    """ 
    plt.scatter(q, sigma_q/q)
    plt.axhline(threshold, color='red', linestyle='--', label="1% Error")
    plt.yscale('log')
    plt.xlabel("Charge")
    plt.ylabel("Fractional Error")
    plt.title("Fractional Error in Charge vs. Charge")
    plt.show()
    
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

    # Use pd.cut to assign each charge value to a bin
    bin_indices = pd.cut(q, bins=bin_edges, labels=False, right=False) #right=False to make bin edges inclusive
    
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

def linear_fit(bin_means):
    """
    Performs a linear regression on the binned charge data.

    Args:
        bin_means (DataFrame): DataFrame with columns ['bin', 'mean_charge', 'sigma_mean']

    Returns:
        tuple: (slope, intercept, slope_error, intercept_error)
    """
    x_data = bin_means["bin"]  # Bin indices (or centers if preferred)
    y_data = bin_means["mean_charge"]
    sigma_y = bin_means["sigma_mean"]  # Use uncertainties as weights

    # Define linear model
    def linear(x, m, b):
        return m * x + b

    # Perform weighted least squares fit
    popt, pcov = curve_fit(linear, x_data, y_data, sigma=sigma_y, absolute_sigma=True)

    # Extract parameters and uncertainties
    slope, intercept = popt
    slope_error, intercept_error = np.sqrt(np.diag(pcov))

    return slope, intercept, slope_error, intercept_error


def plot_binned_data(bin_means, slope, intercept):
    """
    Plots the binned charge data along with the linear regression fit.

    Args:
        bin_means (DataFrame): DataFrame with columns ['bin', 'mean_charge', 'sigma_mean']
        slope (float): Slope of the linear fit
        intercept (float): Intercept of the linear fit
    """
    plt.errorbar(bin_means["bin"], bin_means["mean_charge"], yerr=bin_means["sigma_mean"], fmt='o', mfc='none',label="Binned Data")
    x_fit = np.linspace(min(bin_means["bin"]), max(bin_means["bin"]), 100)
    plt.plot(x_fit, slope * x_fit + intercept, label="Linear Fit", linestyle="--")
    
    plt.xlabel("Bin Index")
    plt.ylabel("Mean Charge")
    plt.legend()
    plt.title("Linear Fit to Binned Charge Data")
    plt.show()
