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
    r = Pos(x)
    theta = in_tan(r)
    return (theta)

def Pos(x, s):
    """
    Computes the position relative to a shift and propagates uncertainty.

    Parameters:
        x, s: Tuples of the form (value, uncertainty).
        x: position we measure
        s: "center" of the ruler we measured on
    
    Returns:
        A tuple (pos, sigma_pos) where:
        - pos is the computed position.
        - sigma_pos is the propagated uncertainty.
    """
    # Extract values and uncertainties
    x_val, sigma_x = x
    s_val, sigma_s = s

    # Compute the function value
    pos = x_val - s_val

    # Compute the propagated uncertainty
    sigma_pos = np.sqrt(sigma_x**2 + sigma_s**2)

    return (pos, sigma_pos)


def Theta(r, l=(6.6, 0.002)):
    """
    Computes Theta as ArcTan(r/l)/2 and propagates uncertainty.
    Note: gives *half* the angle on the ruler, since the law of reflection doubles the actual angle of the beam

    Parameters:
        r, l: Tuples of the form (value, uncertainty).
        r: the position relative to the center of the ruler
        l: distance from device to the wall
    Returns:
        A tuple (theta, sigma_theta) where:
        - theta is the computed angle.
        - sigma_theta is the propagated uncertainty.
    """
    # Extract values and uncertainties
    r_val, sigma_r = r
    l_val, sigma_l = l

    # Compute the function value
    theta = 0.5 * np.arctan(r_val / l_val)

    # Compute the propagated uncertainty
    denom = 4 * (1 + (r_val / l_val))**2 * l_val**2
    sigma_theta = np.sqrt((sigma_r**2 / denom) + ((r_val**2 * sigma_l**2) / (denom * l_val**2)))

    return (theta, sigma_theta)


def I_beam(m, l, w):
    """
    Computes the moment of inertia for an I-beam and propagates uncertainty.

    Parameters:
        m, l, w: Tuples of the form (value, uncertainty).
        
        m (tuple): mass of the beam
        l (tuple): length of beam
        w (tuple): width of beam
    Returns:
        A tuple (I, sigma_I) where:
        - I is the computed moment of inertia.
        - sigma_I is the propagated uncertainty.
    """
    # Extract values and uncertainties
    m_val, sigma_m = m
    l_val, sigma_l = l
    w_val, sigma_w = w

    # Compute the function value
    I = m_val * (1/12) * (l_val**2 + w_val**2)

    # Compute the propagated uncertainty
    term1 = (1/144) * (l_val**2 + w_val**2)**2 * sigma_m**2
    term2 = (1/36) * m_val**2 * l_val**2 * sigma_l**2
    term3 = (1/36) * m_val**2 * w_val**2 * sigma_w**2

    sigma_I = np.sqrt(term1 + term2 + term3)

    return (I, sigma_I)

def ISpheres(m, dss, d):
    """
    Computes the moment of inertia for two spheres and propagates uncertainty.

    Parameters:
        m, dss, d: Tuples of the form (value, uncertainty).
    
    Returns:
        A tuple (I, sigma_I) where:
        - I is the computed moment of inertia.
        - sigma_I is the propagated uncertainty.
    """
    # Extract values and uncertainties
    m_val, sigma_m = m
    dss_val, sigma_dss = dss
    d_val, sigma_d = d

    # Compute the function value
    I = 2 * m_val * ((2/5) * (dss_val / 2)**2 + d_val**2)

    # Compute the propagated uncertainty
    term1 = 4 * ((dss_val**2 / 10) + d_val**2)**2 * sigma_m**2
    term2 = (4/25) * m_val**2 * dss_val**2 * sigma_dss**2
    term3 = 16 * m_val**2 * d_val**2 * sigma_d**2

    sigma_I = np.sqrt(term1 + term2 + term3)

    return (I, sigma_I)

def Itotal(Ib, Is):
    """
    Computes the total moment of inertia and propagates uncertainty.

    Parameters:
        Ib, Is: Tuples of the form (value, uncertainty).
    
    Returns:
        A tuple (I_total, sigma_I_total) where:
        - I_total is the computed total moment of inertia.
        - sigma_I_total is the propagated uncertainty.
    """
    # Extract values and uncertainties
    Ib_val, sigma_Ib = Ib
    Is_val, sigma_Is = Is

    # Compute the function value
    I_total = Ib_val + Is_val

    # Compute the propagated uncertainty
    sigma_I_total = np.sqrt(sigma_Ib**2 + sigma_Is**2)

    return (I_total, sigma_I_total)