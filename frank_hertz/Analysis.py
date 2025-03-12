"""
Script containing functions for analysis of the Frank-Hertz experiment. Heavily lifted from the millikan analysis script.
"""
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
def linear_fit(df):
    """
    Performs a linear regression on the points taken from frank hertz images.

    Args:
        df (DataFrame): DataFrame including columns ['n', 'V', 'sigma_V']

    Returns:
        tuple: (slope, intercept, slope_error, intercept_error)
    """
    x_data = df["n"]  # n indices (or centers if preferred)
    y_data = df["V"]
    sigma_y = df['sigma_V']  # Use uncertainties as weights
    sigma_y = df['big_sigma_v']  # Use uncertainties as weights

    # Define linear model
    def linear(x, m, b):
        return m * x + b

    # Perform weighted least squares fit
    popt, pcov = curve_fit(linear, x_data, y_data, sigma=sigma_y, absolute_sigma=True)

    # Extract parameters and uncertainties
    slope, intercept = popt
    slope_error, intercept_error = np.sqrt(np.diag(pcov))

    return slope, intercept, slope_error, intercept_error


def plot_data(df, slope, intercept):
    """
    Plots the data along with the linear regression fit.

    Args:
        df(DataFrame): DataFrame with columns ['n', 'V', 'sigma_V']
        slope (float): Slope of the linear fit
        intercept (float): Intercept of the linear fit
    """
    plt.errorbar(df["n"], df["V"], yerr=df["sigma_V"], fmt='o', mfc='none',)
    x_fit = np.linspace(min(df["n"]), max(df["n"]), 100)
    plt.plot(x_fit, slope * x_fit + intercept, label="Linear Fit", linestyle="--")
    
    plt.xlabel("n Index")
    plt.ylabel("Voltage")
    plt.legend()
    plt.title("Linear Fit to Frank-Hertz Data")
    plt.show()
