import pandas as pd
import math
import numpy as np
from cavendish.utils.Functions import *
from compton.Functions import *
from millikan.functions import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot_raw(df, angle, x_label="Energy", y_label="Counts"):
    #plots from the dataframe
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    plt.plot(energy, counts, label=f"A{angle}")
    plt.title(f"A{angle} Raw Data")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.show()
    
    
def plot_all(dfs, mins,):
    #plots all the data from the dataframe
    
    for i, angle in enumerate(range(0, 140, 10)):
        plt.axvline(x=mins[i], color='r', linestyle='--')
        plot_raw(dfs[i], angle)
        
def plot_all_S1_fits(dfs):
    #plots all the fits from the dataframe
    #used to also return the peaks, but that was removed to prevent confusion. data returning functionality is provided in run_all_S1_fits
    for i, angle in enumerate(range(0, 140, 10)):
        plot_guassian_fit(dfs[i], i, angle)
    
def run_single_S1_fit(df, i=0, angle=0):
    #essentially a function wrapper for running one at a time
    dfs = [df]
    return run_all_S1_fits(dfs)

def run_all_S1_fits(dfs):
    #runs all the fits from the dataframe, returns dataframe with statistics for the peaks
    rows= []
    for i, df in enumerate(dfs):
        popt, pcov = guassian_fit(df, i, i*10)
        mean, sigma = popt[1], popt[2]
        Unc_mean, Unc_sigma = np.sqrt(np.diag(pcov)[1:3])
        rows.append({"Angle": i*10, "Mean": mean, "Sigma": sigma, "Unc Mean": Unc_mean, "Unc Sigma": Unc_sigma})
    peaks = pd.DataFrame(rows)
    return peaks
def trim_S1_dfs(dfs, mins, maxs):
    #trims the dataframes to the minimum energy
    trimmed = []
    for i, df in enumerate(dfs):
        trimmed.append(trim_df(df, mins[i], maxs[i]))
    return trimmed

def trim_df(df, min_energy, max_energy=None):
    #trims the dataframe to the minimum energy
    if max_energy is None:
        return df[df["Energy (keV)"] > min_energy]
    else:
        return df[(df["Energy (keV)"] > min_energy) & (df["Energy (keV)"] < max_energy)]

def guassian_fit(df, i=0, angle=0, p0_overide = None):
    #fits a guassian to the data, returns parameters and cov matrix
    #should figure out why I added i to this, probably for compatability with an enumeration function
    if type(angle) != int:
        angle = 0
        print("Angle not int, defaulting to 0")
    mins = minimum_energy_S1()
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    if p0_overide != None:
        p0 = p0_overide
    else:
        p0=[safe_divide(30000, angle), mins[i]+ 50, 100]
    popt, pcov = curve_fit(gaussian, energy, counts, p0=p0)
    return popt, pcov

def plot_guassian_fit(df, i=0, angle=0, p0_overide = None):
    #plots the guassian fit
    #used to return data, now does not, as to prevent confusion. That was moved to run_all_S1_fits
    popt, pcov= guassian_fit(df, i, angle, p0_overide=p0_overide)
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    plt.plot(energy, counts, label=f"A{angle}")
    plt.plot(energy, gaussian(energy, *popt), label=f"A{angle} Fit")
    plt.title(f"A{angle} Guassian Fit")
    plt.xlabel("Energy")
    plt.ylabel("Counts")
    plt.legend()
    plt.show()

def error_fit(df):
    #fits an error function to a dataframe containing the ROI for the compton edge
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    popt, pcov = curve_fit(fit_erf, energy, counts, p0=[max(counts), min(counts), 1, 475])
    return popt, pcov

def plot_error_fit(df):
    #plots the error fit
    popt, pcov = error_fit(df)
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    plt.plot(energy, counts, label="Raw Data")
    plt.plot(energy, fit_erf(energy, *popt), label="Error Fit")
    plt.title("Error Fit to Data")
    plt.xlabel("Energy")
    plt.ylabel("Counts")
    plt.legend()
    plt.show()

def logistic_fit(df):
    #fits a logistic function to the data, returns parameters and cov matrix
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    popt, pcov = curve_fit(logistic, energy, counts, p0=[max(counts), min(counts), 1, 475])
    return popt, pcov

def plot_logistic_fit(df):
    #plots the logistic fit
    popt, pcov = logistic_fit(df)
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    plt.plot(energy, counts, label="Raw Data")
    plt.plot(energy, logistic(energy, *popt), label="Logistic Fit")
    plt.title("Logistic Fit to Data")
    plt.xlabel("Energy")
    plt.ylabel("Counts")
    plt.legend()
    plt.show()

def linear_fit_plot(peaks):
    """creates a fit from an expanded peaks dataframe, plots the fit and reports parameters with uncertainty
    """
    x = peaks["1-cos(theta)"]
    y = peaks["E/E'"]
    popt, pcov = curve_fit(linear, x, y, p0=[1, 0])
    slope, unc_slope =popt[0], np.sqrt(np.diag(pcov)[0])
    intercept, unc_intercept = popt[1], np.sqrt(np.diag(pcov)[1])
    plt.errorbar(x, y, yerr=peaks["Unc E/E'"], fmt='o', mfc='none')
    x_fit = np.linspace(min(x), max(x), 100)
    plt.plot(x_fit, linear(x_fit, slope, intercept), label="Linear Fit", linestyle="--")
    plt.xlabel("1-cos(theta)")
    plt.ylabel("E/E'")
    plt.legend()
    plt.title("Linear Fit to Expanded Peaks Data")
    plt.show()
    return slope, unc_slope, intercept, unc_intercept



def minimum_energy_positron_peaks():
    return [185, 240, 460]

def minimum_energy_S1():
    #returns handmade list of reasonable energy cutoffs to enable fits. note: only for S1 runs
    return [650, 650, 600, 500, 450, 400, 350, 300, 280, 250, 240, 220, 220, 200]

def maximum_energy_S1():
    #returns handmade list of reasonable energy cutoffs to enable fits. note: only for S1 runs
    return [710, 710, 700, 620, 575, 520, 475, 420, 380, 340, 320, 280, 270, 250]