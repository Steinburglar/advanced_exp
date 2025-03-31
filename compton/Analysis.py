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
    rows= []
    for i, angle in enumerate(range(0, 140, 10)):
        mean, sigma = plot_guassian_fit(dfs[i], i, angle)
        rows.append({"Angle": angle, "Mean": mean, "Sigma": sigma})
    peaks = pd.DataFrame(rows)
    
def run_all_S1_fits(dfs):
    #runs all the fits from the dataframe
    rows= []
    for i, df in enumerate(dfs):
        _, mean, sigma = guassian_fit(df, i, i*10)
        rows.append({"Angle": i*10, "Mean": mean, "Sigma": sigma})
    peaks = pd.DataFrame(rows)
    return peaks
def trim_S1_dfs(dfs, mins):
    #trims the dataframes to the minimum energy
    trimmed = []
    for i, df in enumerate(dfs):
        trimmed.append(df[df["Energy (keV)"] > mins[i]])
    return trimmed

def guassian_fit(df, i=None, angle=None):
    #fits a guassian to the data
    mins = minimum_energy_S1()
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    popt, pcov = curve_fit(gaussian, energy, counts, p0=[safe_divide(30000, angle), mins[i]+ 50, 100])
    return popt

def plot_guassian_fit(df, i=None, angle=None):
    #plots the guassian fit, returns the mean and sigma
    popt = guassian_fit(df, i, angle)
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    plt.plot(energy, counts, label=f"A{angle}")
    plt.plot(energy, gaussian(energy, *popt), label=f"A{angle} Fit")
    plt.title(f"A{angle} Guassian Fit")
    plt.xlabel("Energy")
    plt.ylabel("Counts")
    plt.legend()
    plt.show()
    return popt[1], popt[2]

def minimum_energy_S1():
    #returns handmade list of reasonable energy cutoffs to enable fits. note: only for S1 runs
    return [550, 550, 520, 450, 400, 350,320, 300, 280, 250, 220,210, 200, 190]