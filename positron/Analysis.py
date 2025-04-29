import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from compton.Analysis import *
from compton.Dataloader import *
from compton.Functions import *
from positron.Functions import *
from cavendish.utils.Functions import *
from compton.Functions import *
from millikan.functions import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def plot(df, title, x_label="Energy", y_label="Counts"):
    #plots from the dataframe
    energy = df["Energy (keV)"].to_numpy()
    counts = df["Counts"].to_numpy()
    plt.plot(energy, counts,)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.show()


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
        #find the normalization area for overlap
        x = np.arange(-6, 7, 1)
        discrete_overlap = [area_overlap(angle) for angle in x]
        area = np.sum(discrete_overlap)
        #normalize the overlap
        
        x = np.linspace(-6, 6, 100)
        normalized_overlap = [area_overlap(angle)/area for angle in x]
        plt.plot(x, normalized_overlap, label="Overlap", color="red")
        plt.legend()
    
    plt.show()
    
    

def normalize_counts(df):
    """Normalizes the counts in a DataFrame to the area.

    Args:
        df (pd.DataFrame): dataframe with columns "angle" and "counts".

    Returns:
        pd.DataFrame: input dataframe with new collumn of normalized counts, and normalized errors
    """
    area = df["Counts"].sum()
    df["normalized_counts"] = df["Counts"] / area
    df["normalized_error"] = df["Error"] / area
    return df


def total_counts(dfs):
    #takes in a dictionary of dataframes for each angle, returns a dataframe with the total counts and errors at each angle
    counts = {angle: df["Counts"].sum() for angle, df in dfs.items()} #create a dictionary from the dfs
    counts = pd.DataFrame.from_dict(counts, orient='index', columns=["Counts"]) # create a dataframe from the dictionary
    counts.index.name = "Angle" # set the index name to "Angle"

    counts["Error"] = np.sqrt(counts["Counts"])
    counts = normalize_counts(counts)
    return counts

def total_roi_counts(dfs):
    #takes in a dictionary of dataframes for each ROI, returns a dataframe with the total counts at each ROI
    
    
    
    counts = {roi: df["Counts"].sum() for roi, df in dfs.items()} #create a dictionary from the dfs
    counts = pd.DataFrame.from_dict(counts, orient='index', columns=["Counts"]) # create a dataframe from the dictionary
    counts.index.name = "ROI" # set the index name to "Angle"

    counts["Error"] = np.sqrt(counts["Counts"])
    return counts

def run_all_peak_fits(dfs, p0_overide=None):
    #runs fits for all ROI peaks from Atotal dataframe, returns dataframe with statistics for the peaks
    #args: dfs. this is a DICTIONARY of trimmed dataframes for each  ROI of the form {"ROI": dataframe}
    rows= []
    for roi, df in dfs.items():
        p0_overide = [1000, 250, 20] if roi == "Annihilation" else p0_overide
        p0_overide = [500, 260, 20] if roi == "Cesium" else p0_overide
        p0_overide = [175, 490, 20] if roi == "High" else p0_overide
        popt, pcov = guassian_fit(df, p0_overide=p0_overide)
        mean, sigma = popt[1], popt[2]
        Unc_mean, Unc_sigma = np.sqrt(np.diag(pcov)[1:3])
        rows.append({"ROI": roi, "Mean": mean, "Sigma": sigma, "Unc Mean": Unc_mean, "Unc Sigma": Unc_sigma})
    peaks = pd.DataFrame(rows)
    peaks.set_index("ROI", inplace=True)
    return peaks

def fit_plot_calibration(df, knowns):
    #fits a line to the three known energy peaks to rescale the energy axis.
    x = df["Mean"].to_numpy()
    y = knowns
    popt, pcov = curve_fit(linear, x, y,)
    slope, intercept = popt
    Unc_slope, Unc_intercept = np.sqrt(np.diag(pcov))
    plt.plot(x, y, "o", label="Data")
    x_ = np.linspace(min(x), max(x), 100)
    plt.plot(x_, linear(x_, slope, intercept), label="Fit")
    plt.show
    return slope, intercept, Unc_slope, Unc_intercept
    
    
def recalibrate_energy(df, mins, maxs):
    #identifies peaks for known energy peaks, performs fit and recalibrates energy.
    #returns dataframe with recalibrated energy
    #args: df, dataframe with columns "Energy (keV)" and "Counts"
    #       mins, list of minimum energies for each ROI
    #       maxs, list of maximum energies for each ROI
    
    dfs = [trim_df(df, min_, max_) for min_, max_ in zip(mins, maxs)]
    if len(dfs) == 2:
        ROIs = {"Annihilation": dfs[0], "High": dfs[1]}
        peaks = run_all_peak_fits(ROIs, p0_overide=[100, 250, 20])
        display(peaks)
        slope, intercept, sig_slope, sig_int = fit_plot_calibration(peaks, [511, 1274])
    elif len(dfs) == 3:
        ROIs = {"Annihilation": dfs[0], "Cesium": dfs[1], "High": dfs[2]}
        peaks = run_all_peak_fits(ROIs)
        slope, intercept, sig_slope, sig_int = fit_plot_calibration(peaks, [511, 662, 1274])
    else:
        raise ValueError("Invalid number of ROIs")
    print(slope, intercept)
    print(sig_slope, sig_int)
    df["Energy (keV)"] = linear(df["Energy (keV)"], slope, intercept) #recalibrate energy
    return df

def time_normalize(rates):
    #takes dataframes with rates at each ROI, normalizes by time of observation. should only be used once in analysis, in part 2 where we asses window inter val
    rates["n1"] = rates["n1"]*1200 /(6757.0260)
    rates["n1_error"] = rates["n1_error"]*1200 /(6757.0260)
    rates["n2"] = rates["n2"]*1200 /(300)
    rates["n2_error"] = rates["n2_error"]*1200 /(300)
    rates["N_acc"] = rates["N_acc"]
    rates["N_acc_error"] = rates["N_acc_error"]
    return rates

def window_error(n1, n2, n1_err, n2_err, N_acc, N_accerr):
    #calculates the error in the window interval
    #n1, n2, n1err, n2err, N_acc, N_accerr are all scalars
    #returns the error in time window
    term1 = np.square(N_acc*n1_err/(n1*n1*n2))
    term2 = np.square(N_acc*n2_err/(n1*n2*n2))
    term3 = np.square(N_accerr/(n1*n2))
    return np.sqrt(term1 + term2 + term3)
