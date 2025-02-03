import numpy as np
import numpy as np
import math

"""
script containing necessary functions for analysis, once the Paramters from fitting the data have been obtained, as well as the function to be fit itself
With the exception of Damped Oscillation, all arguments should be duples of (value, Uncertainty). this allows us to paralellize our uncertainty calculations.
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
    
    
def con_rad(x):
    """

    Args:
        x (duple): (position measurement, uncertainty on that position)

    Returns:
        _type_: (theta measurement, uncertianty on that theta measurement)
    """\
    #mathematical composed function to convert from position to radians
    #takes tuples of (value, uncertainty) 
    r = sub_s(x)
    theta = in_tan(r)
    return (theta)

def sub_s(x,s=(1, 0.002)):
    #takes tuples of (value, uncertainty) 
    return (x[0]-s[0], x[1] + s[1])

def in_tan(r, d=(6.6, 0.002)): #needs fixing once we get uncertainty expression
    #takes tuples of (value, uncertainty) 
    #note: also halves the angle, in order to give the angle of rotation the the pendulum, rather than the angle reflected. 
    theta = math.atan(r[0]/d[0])/2
    top1 = (r[0]**2)*(d[1]**2)
    top2 = (d[0]**2)*(r[1]**2)
    bot = 4*(r[0]**2 + d[0]**2)**2
    sig_theta = np.sqrt((top1+top2)/bot)
    return (theta, sig_theta)


def I_beam(mb, lb, wb):
    """calculates Moment of inertia for the beam

    Args:
        mb (tuple): mass of the beam
        lb (tuple): length of beam
        wb (tuple): width of beam
    """