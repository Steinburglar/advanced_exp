"""
Lucas Steinberger
1/31/2025
Python script containing  functions to for loading and converting data from Casey, Lucas, and Cooke's Cavendish lab in PHYS 064.
Goal is to provide utility to load csv and return 4 arrays, each containing the data points with corresponding uncertainties for a measurement,
with position converted to radians of beam rotation.
"""
import pandas as pd
import math


def load_convert(path, sig_x = 0.002, num = 4 ):
    """ Main function for dataloader.
    position is converted to radians of beam rotation. Note that that means the values returned are *half* of the angle produced by the laser, since the law of reflection
    will double the angle between the laser rays.

    Args:
        path (_str_): Path to the location of the .csv file containing measurements of time and position (x) in cm. 
        

    Returns:
        _list_: list of data arrays:  [data_aray1, data_aray2, data_array3 ....]
            
            _data_array_: numoy array containing 3 collumns:
                Time    | angle  | uncertainty in angle
        
    """
    df = pd.read_csv(path, header=[0,1])
    dfs = clean_split(df, num, sig_x)
    convert_rad(dfs)
    return dfs


    
    


def clean_split(df, num, sig_x):
    """takes in a dataframe, cleans it to keep only the Time (sec) and Position (cm) columns, and splits into 4 dataframes, one for each measurement. 
    arguments:
        df: dataframe to be cleaned and split
        num: number of groups
    returns:
        List of dataframes, one for each measurement. 
    """
    df = df.loc[:, pd.IndexSlice[:, ["Time (sec)", "Position (cm)"]]]
    dfs = [df.loc[:, pd.IndexSlice[f"Measurement {i}", :]].copy() for i in range(1, num+1)]
    del df #for memory
    for i, _df in enumerate(dfs):
        _df.loc[:, (f"Measurement {i+1}", "Uncertainty (m)")] = sig_x
        print(_df.columns)
        _df.rename(columns={"Position (cm)": "Position (m)"}, level = 1, inplace = True)
        _df.loc[:, _df.columns.get_level_values(1) == 'Position (m)'] /= 100
        

    
    
    return dfs
    
def convert_rad(dfs):
    for i, df in enumerate(dfs):
        measure = df.columns.get_level_values(0)[0]
        df[[(measure, "Radians"), (measure, "Uncertainty (rad)")]] = \
            df.apply(lambda row: pd.Series(
                con_rad((row[(measure, "Position (m)")],row[(measure, "Uncertainty (m)")])),            \
                    index=[(measure, "Radians"), (measure, "Uncertainty (rad)")]), axis = 1,)
        df.drop([(measure, "Position (m)"), (measure, "Uncertainty (m)")], axis = 1, inplace=True)
    return

def con_rad(x):
    r = sub_s(x)
    theta = in_tan(r)
    return (theta)


    
def sub_s(x,s=(1, 0.002)):
    return (x[0]-s[0], x[1] + s[1])

def in_tan(r, d=(8.5, 0.002)):
    #note: also halves the angle, in order to give the angle of rotation the the pendulum, rather than the angle reflected. 
    theta = math.atan(r[0]/d[0])/2
    sig_theta = r[1]/2
    return (theta, sig_theta)