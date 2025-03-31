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
    
def run_all_S1_fits(dfs):
    #runs all the fits from the dataframe
    rows= []
    for i, df in enumerate(dfs):
        popt, pcov = guassian_fit(df, i, i*10)
        mean, sigma = popt[1], popt[2]
        Unc_mean, Unc_sigma = np.diag(pcov)[1:3]
        rows.append({"Angle": i*10, "Mean": mean, "Sigma": sigma, "Unc Mean": Unc_mean, "Unc Sigma": Unc_sigma})
    peaks = pd.DataFrame(rows)
    return peaks
def trim_S1_dfs(dfs, mins):
    #trims the dataframes to the minimum energy
    trimmed = []
    for i, df in enumerate(dfs):
        trimmed.append(df[df["Energy (keV)"] > mins[i]])
    return trimmed

def trim_df(df, min_energy, max_energy=None):
    #trims the dataframe to the minimum energy
    if max_energy is None:
        return df[df["Energy (keV)"] > min_energy]
    else:
        return df[(df["Energy (keV)"] > min_energy) & (df["Energy (keV)"] < max_energy)]

def guassian_fit(df, i=0, angle=0):
    #fits a guassian to the data, returns parameters and cov matrix
    mins = minimum_energy_S1()
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    popt, pcov = curve_fit(gaussian, energy, counts, p0=[safe_divide(30000, angle), mins[i]+ 50, 100])
    return popt, pcov

def plot_guassian_fit(df, i=None, angle=None):
    #plots the guassian fit
    #used to return data, now does not, as to prevent confusion. That was moved to run_all_S1_fits
    popt, pcov= guassian_fit(df, i, angle)
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    plt.plot(energy, counts, label=f"A{angle}")
    plt.plot(energy, gaussian(energy, *popt), label=f"A{angle} Fit")
    plt.title(f"A{angle} Guassian Fit")
    plt.xlabel("Energy")
    plt.ylabel("Counts")
    plt.legend()
    plt.show()

def minimum_energy_S1():
    #returns handmade list of reasonable energy cutoffs to enable fits. note: only for S1 runs
    return [550, 550, 520, 450, 400, 350,320, 300, 280, 250, 220,210, 200, 190]