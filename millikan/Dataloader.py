"""
Lucas Steinberger
1/31/2025
Python script containing  functions to for loading and converting data from Casey, Lucas, and Cooke's Millikan lab in PHYS 064.
"""
import pandas as pd
import math
import numpy as np
from cavendish.utils.Functions import *
from millikan.functions import *

def load_velocity_data(csv_path):
    """
    Loads velocity data from a CSV file into a nested dictionary.
    
    Args:
        csv_path (str): Path to the CSV file.
        
    Returns:
        dict: Nested dictionary with BBL as keys, Drop numbers as subkeys, and lists of velocities as values.
    """
    df = pd.read_csv(csv_path, header=[0, 1])  # Read with two header rows
    
    velocity_data = {}
    
    # Extract unique BBL labels from the first header row
    unique_bbls = df.columns.get_level_values(0).unique()
    
    for bbl in unique_bbls:
        velocity_data[bbl] = {}
        
        # Filter columns belonging to this BBL
        bbl_columns = df[bbl]
        
        for drop in bbl_columns.columns:  # Second-level headers
            velocity_data[bbl][drop] = bbl_columns[drop].dropna().to_numpy()*(10**-2)# Store non-null values, cm to meters correction included here

    return velocity_data

def load_RV_data(csv_path):
    """
    Loads velocity data from a CSV file into a nested dictionary.
    
    Args:
        csv_path (str): Path to the CSV file.
        
    Returns:
        dict: Nested dictionary with BBL as keys, Drop numbers as subkeys, and lists of velocities as values.
    """
    df = pd.read_csv(csv_path, header=[0])  # Read with two header rows

    

    return df

def main_dataframe(vel_path, volt_path):
    """main function of data loader. Loads and cleans velocities with load_velocity_data. then uses these velocities to populate a 2d dataframe with entries for each droplet.
    
    

    Args:
        vel_path (): _description_
        volt_path (_type_): _description_
    """

    bbl_dict = load_velocity_data(vel_path)
    df = init_df(bbl_dict, volt_path)
    df = expand_df(df)
    return df
    


def init_df(bbl_dict, volt_path):
    """initializes and builds the first phase of the dataframe containing droplets as rows, and relevant data as collumns. does NOT do any calculations of charge or temperature. 

    Args:
        bbl_dict (_dict_): nested dictionary of BBl's and their drops. Should be taken as the output of load_velocity_data. Note: drops are ordered cumulatively,
        so the first drop in BBL 2 is named "drop 8", and so on
        volt_path (_type_): path the the csv containing voltage and resistance measurement for each BBL

    Returns:
        _DataFrame_: and Pandas dataframe with one row per droplet, named "drop N", and columns for the initial data. 
    """
    df = pd.DataFrame(columns=["v_rise", "sigma_v_rise", "v_fall", "sigma_v_fall", "Volts", "sigma_Volts", "resistance"])
    volt_df = pd.read_csv(volt_path)
    for bbl, drops in bbl_dict.items():
        bblV = bbl + "V"
        bblR = bbl + "R"
        V_low = volt_df[bblV][0]
        V_hi = volt_df[bblV][1]
        V_mean = np.mean([V_low, V_hi])
        sigV = (V_hi-V_low)/2
        Res = np.mean([volt_df[bblR][0], volt_df[bblR][1]])
        for drop, velocities in drops.items():
            (v_rise, sigma_v_rise), (v_fall, sigma_v_fall) = process_droplet_velocities(velocities)
            df.loc[drop] = [v_rise, sigma_v_rise, v_fall, sigma_v_fall, V_mean, sigV, Res]     
    df = df.dropna()
    
    return df

def calculate_charge(row):
        """small function put together functions that serve to calculate A. this is designed to be used by df.apply()"""
        a_ = a_a((row["viscosity"], row["sigma_viscosity"]))
        q, sigma_q = Q((row["v_fall"], row["sigma_v_fall"]), (row["v_rise"], row["sigma_v_rise"]), (row["Efield"], row["sigma_Efield"]), (row["pressure"], 0), (row["viscosity"], row["sigma_viscosity"]), 8.2*(10**-3))
        return (q, sigma_q)
    
    #(Efield, Vfall, Vrise, eta, rho, g, b, pressure):
    
def calculate_charge_2(row):
    q, sigma_q = drop_charge_and_uncertainty((row["Efield"], row["sigma_Efield"]), (row["v_fall"], row["sigma_v_fall"]), (row["v_rise"], row["sigma_v_rise"]), row["viscosity"], row["pressure"])
    return (q, sigma_q)

def expand_df(df):
    """function to fill in the calculated and extra nformation for each drop in the main dataframe. 

    Args:
        df (_DataFrame_): dataframe containing a row for each drop, and initialized with some already existing columns. should be made by init_df
    """

    df["Temperature"]= df.apply(lambda row: resistance_to_temperature((row["resistance"])), axis = 1, result_type="expand")
    df["pressure"]= inches_to_pascals(30.28)
    df[["spacing", "sigma_spacing"]]= [0.00768167, 0.00000752773] #need to add actual spacing uncertainty
    df[["Efield", "sigma_Efield"]] = df.apply(
        lambda row: efield((row["Volts"], row["sigma_Volts"]), (row["spacing"], row["sigma_spacing"])), 
        axis=1, result_type="expand"
    )
    df[["viscosity", "sigma_viscosity"]] = df.apply(lambda row: viscosity((row["Temperature"], 0)), axis = 1, result_type="expand")

    df[["q", "sigma_q"]] = df.apply(calculate_charge, axis = 1, result_type = "expand")
    
    return df

def load_class_charges(class_path):
    """loads class data as a dataframe, extracts charges and uncertainty in charges. 

    Args:
        class_path (str): path to class csv from sliwa
    """
    # Load only the first two columns
    df = pd.read_csv(class_path, header=0)
    df = df[df['delta_t_fall[s]'] > 0]
    
    df[df.select_dtypes(include=['number']).columns] = df.select_dtypes(include=['number']).abs()
    df = df.rename(columns={"temperature[K]": "Temperature", "air_pressure[Pa]": "pressure"})
    df[["v_fall", "sigma_v_fall"]]= df.apply(lambda row: velocity((row["delta_x_fall[m]"], row["sigma_delta_x_fall[m]"]), row["delta_t_fall[s]"],), axis = 1, result_type="expand")
    df[["v_rise", "sigma_v_rise"]]= df.apply(lambda row: velocity((row["delta_x_rise[m]"], row["sigma_delta_x_rise[m]"]), row["delta_t_rise[s]"],), axis = 1, result_type="expand")
    #df[["spacing", "sigma_spacing"]]= [0.00768167, 0.00000752773] #need to add actual spacing uncertainty
    df[["spacing", "sigma_spacing"]]= [0.00767, 0.00002]
    df[["Efield", "sigma_Efield"]]= df.apply(lambda row: efield((row["voltage[V]"], row["sigma_voltage[V]"]), (row["spacing"], row["sigma_spacing"])), axis = 1, result_type="expand")
    df[["viscosity", "sigma_viscosity"]] = df.apply(lambda row: viscosity((row["Temperature"], 0)), axis = 1, result_type="expand")
    df[["q", "sigma_q"]] = df.apply(calculate_charge_2, axis = 1, result_type = "expand")
    return df