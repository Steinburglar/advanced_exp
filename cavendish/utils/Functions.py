import numpy as np

"""
script containing necessary functions for analysis, once the paramters from fitting the data have been obtained, as well as the function to be fit itself
"""

def damped_oscillation(time, theta0, amp, period, delta, b):
    """
    Function to be fit to the data. operates on a point-wise basis.
    
    Arguments:
    time (_float_): time value
    theta0 (_float_): equilibrium value of theta
    amp (_float_): amplitude of osscilations
    period (_float_):
    delta (_float_): 
    b (_float_)
    """
    return theta0 + amp*np.sin((2(np.pi)/period)*time+delta)*np.exp(-b*time)
    
