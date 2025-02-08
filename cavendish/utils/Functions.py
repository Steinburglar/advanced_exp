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
    amp (_float_): amplitude of oscilations
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
    theta = Theta(r)
    return (theta)

def Pos(x, s=(0.86917, 0.0012)):
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


def Theta(r, l=(6.6385, 0.021)):
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
    denom = (1 + (r_val / l_val))**2 * l_val**2
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
        Ib
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

def Lambda(i, T):
    """
    Computes Lambda as 4 * (pi^2) * i / (T^2) and propagates uncertainty.

    Parameters:
        i, T: Tuples of the form (value, uncertainty).
        i: total moment of inertia
        T: period of oscillation (seconds)
    
    Returns:
        A tuple (Lambda_value, sigma_Lambda) where:
        - Lambda_value is the computed function value.
        - sigma_Lambda is the propagated uncertainty.
    """
    # Extract values and uncertainties
    i_val, sigma_i = i
    T_val, sigma_T = T

    # Compute the function value
    Lambda_value = (4 * np.pi**2 * i_val) / (T_val**2)

    # Compute the propagated uncertainty
    term1 = (16 * np.pi**4 * sigma_i**2) / (T_val**4)
    term2 = (64 * np.pi**4 * i_val**2 * sigma_T**2) / (T_val**6)
    sigma_Lambda = np.sqrt(term1 + term2)

    return (Lambda_value, sigma_Lambda)


def DeltaTheta(theta1, theta2):
    """
    Computes DeltaTheta as (theta1 - theta2) / 2 and propagates uncertainty.

    Parameters:
        theta1, theta2: Tuples of the form (value, uncertainty).
    
    Returns:
        A tuple (DeltaTheta_value, sigma_DeltaTheta) where:
        - DeltaTheta_value is the computed function value.
        - sigma_DeltaTheta is the propagated uncertainty.
    """
    # Extract values and uncertainties
    theta1_val, sigma_theta1 = theta1
    theta2_val, sigma_theta2 = theta2

    # Compute the function value
    DeltaTheta_value = (theta1_val - theta2_val) / 2

    # Compute the propagated uncertainty
    sigma_DeltaTheta = np.sqrt(sigma_theta1**2 / 4 + sigma_theta2**2 / 4)

    return (DeltaTheta_value, sigma_DeltaTheta)

import numpy as np

def DeltaTau(lambda_, deltheta):
    """
    Computes DeltaTau as lambda * deltheta and propagates uncertainty.

    Parameters:
        lambda_, deltheta: Tuples of the form (value, uncertainty).
    
    Returns:
        A tuple (DeltaTau_value, sigma_DeltaTau) where:
        - DeltaTau_value is the computed function value.
        - sigma_DeltaTau is the propagated uncertainty.
    """
    # Extract values and uncertainties
    lambda_val, sigma_lambda = lambda_
    deltheta_val, sigma_deltheta = deltheta

    # Compute the function value
    DeltaTau_value = lambda_val * deltheta_val

    # Compute the propagated uncertainty
    sigma_DeltaTau = np.sqrt(
        (deltheta_val**2) * (sigma_lambda**2) + 
        (lambda_val**2) * (sigma_deltheta**2)
    )

    return (DeltaTau_value, sigma_DeltaTau)


def G_first_order(deltau, b, d, M, m):
    """
    Computes G_first_order and propagates uncertainty.

    Parameters:
        deltau, b, d, M, m: Tuples of the form (value, uncertainty).
        deltau: difference in torques in both positions
        b: distance from the center of lead sphere to center of tungsten sphere (different than the value in the book bc of window change)
        d: distance from the center of the beam to center of the small lead sphere
        M: mass of large tungsten sphere
        m: mass of small lead sphere
    Returns:
        A tuple (G_value, sigma_G) where:
        - G_value is the computed function value.
        - sigma_G is the propagated uncertainty.
    """
    # Extract values and uncertainties
    deltau_val, sigma_deltau = deltau
    b_val, sigma_b = b
    d_val, sigma_d = d
    M_val, sigma_M = M
    m_val, sigma_m = m

    # Compute function value
    G_value = (deltau_val * (b_val**2)) / (4 * m_val * M_val * d_val)

    # Compute uncertainty propagation
    sigma_G = np.sqrt(
        ((b_val**4) * sigma_deltau**2) / (16 * m_val**2 * M_val**2 * d_val**2) +
        ((4 * b_val**2 * deltau_val**2) * sigma_b**2) / (16 * m_val**2 * M_val**2 * d_val**2) +
        ((b_val**4 * deltau_val**2) * sigma_d**2) / (16 * m_val**2 * M_val**2 * d_val**4) +
        ((b_val**4 * deltau_val**2) * sigma_M**2) / (16 * m_val**2 * M_val**4 * d_val**2) +
        ((b_val**4 * deltau_val**2) * sigma_m**2) / (16 * m_val**4 * M_val**2 * d_val**2)
    )

    return (G_value, sigma_G)

def Theta1(b, d):
    """
    Computes Theta1, the angle between the tungsten sphere and the further lead sphere and propagates uncertainty.

    Parameters:
        b, d: Tuples of the form (value, uncertainty).
        b: distance from the center of lead sphere to center of tungsten sphere (different than the value in the book bc of window change)
        d: distance from the center of the beam to center of the small lead sphere
    Returns:
        A tuple (theta1_value, sigma_theta1) where:
        - theta1_value is the computed function value.
        - sigma_theta1 is the propagated uncertainty.
    """
    # Extract values and uncertainties
    b_val, sigma_b = b
    d_val, sigma_d = d

    # Compute function value
    theta1_value = np.arctan(b_val / (2 * d_val))

    # Compute uncertainty propagation
    factor = 1 + (b_val**2) / (4 * d_val**2)
    sigma_theta1 = np.sqrt(
        (sigma_b**2) / (4 * factor**2 * d_val**2) +
        (b_val**2 * sigma_d**2) / (4 * factor**2 * d_val**4)
    )

    return (theta1_value, sigma_theta1)

def GSecond(deltau, b, d, M, m, theta1, r):
    """
    Computes GSecond and propagates uncertainty.

    Parameters:
        deltau, b, d, M, m, theta1, r: Tuples of the form (value, uncertainty).
    
    Returns:
        A tuple (G, sigma_G) where:
        - G is the computed value of the function.
        - sigma_G is the propagated uncertainty.
    """
    # Extract values and uncertainties
    deltau_val, sigma_deltau = deltau
    b_val, sigma_b = b
    d_val, sigma_d = d
    M_val, sigma_M = M
    m_val, sigma_m = m
    theta1_val, sigma_theta1 = theta1
    r_val, sigma_r = r

    # Compute the function value
    factor = 1 / (b_val**2) - np.sin(theta1_val) / r_val**2
    G = (deltau_val / (4 * M_val * m_val * d_val)) * (factor ** -1)

    # Compute the propagated uncertainty
    denom = 16 * d_val**2 * M_val**2 * m_val**2 * factor**2
    term1 = sigma_deltau**2 / denom
    term2 = (deltau_val**2 * sigma_d**2) / (4 * d_val**6 * M_val**2 * m_val**2 * factor**4)
    term3 = (deltau_val**2 * sigma_M**2) / (16 * d_val**2 * M_val**4 * m_val**2 * factor**2)
    term4 = (deltau_val**2 * sigma_m**2) / (16 * d_val**2 * M_val**2 * m_val**4 * factor**2)
    term5 = (deltau_val**2 * sigma_b**2) / (16 * d_val**2 * M_val**2 * m_val**2 * b_val**4 * factor**4)
    term6 = (np.cos(theta1_val)**2 * deltau_val**2 * sigma_theta1**2) / (16 * d_val**2 * M_val**2 * m_val**2 * factor**4 * r_val**4)
    term7 = (np.sin(theta1_val)**2 * deltau_val**2 * sigma_r**2) / (4 * d_val**2 * M_val**2 * m_val**2 * factor**4 * r_val**6)

    sigma_G = np.sqrt(term1 + term2 + term3 + term4 + term5 + term6 + term7)

    return (G, sigma_G)
