"""This script contains functions for use in other parts of the analysis.
"""

import numpy as np
import math

def process_droplet_velocities(velocities):
    """cleans velocity data by fiding average rising and falling velocities, using standard deviation as uncertainty.

    Args:
        velocities (_list_): list of velocities

    3 Example usage:
    velocities = [10, 12, 15, 100, -5, -8, -12, -50, 7, 9, -7, -9, 14, 16, -6, -20]
    rise_stats, fall_stats = process_droplet_velocities(velocities)

    print("Rising droplets (mean velocity, std dev):", rise_stats)
    print("Falling droplets (mean velocity, std dev):", fall_stats)

    Returns:
        _tuple_: tuple of tuples.
    """
    
    
    # Separate rising (positive) and falling (negative) velocities
    v_rise = [v for v in velocities if v > 0]  # Rising droplets
    v_fall = [v for v in velocities if v < 0]  # Falling droplets

    def remove_outliers(data):
        if len(data) < 2:
            return []  # Not enough data to compute IQR
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        return [v for v in data if lower_bound <= v <= upper_bound]

    # Trim outliers from rising and falling velocities
    v_rise_trimmed = remove_outliers(v_rise)
    v_fall_trimmed = remove_outliers(v_fall)

    def compute_mean_stddev(data):
        if len(data) == 0:
            return (None, None)  # Return None if no valid data remains
        return (np.mean(data), 2*np.std(data, ddof=1))  # ddof=1 for sample std dev

    # Compute mean and standard deviation for rising and falling droplets
    rise_stats = compute_mean_stddev(v_rise_trimmed)
    neg_fall_stats = compute_mean_stddev(v_fall_trimmed)
    if rise_stats == (None, None):
        print("Droplet found without enough values")
    fall_stats= np.abs(neg_fall_stats)
    return rise_stats, fall_stats


def a_a(n):
    """
    Computes a (as written in the lab notebook. given n, ro, and g.
    
    Parameters:
        n (tuple): (value, uncertainty) for n
        ro (tuple): (value, uncertainty) for ro
        g (tuple): (value, uncertainty) for g
    
    Returns:
        tuple: (propagated value, propagated uncertainty)
    """
    # Extract values and uncertainties
    n_val, sigma_n = n
    ro_val = 886
    g_val= 9.8

    # Compute function value
    A_val = np.sqrt((9 * n_val**3) / (2 * ro_val * g_val))

    # Compute uncertainty propagation
    sigma_A = (9 / (2 * np.sqrt(2))) * np.sqrt(
        (sigma_n**2 * n_val) / (g_val * ro_val)
    )

    return A_val, sigma_A


def Q(vf, vr, Efield, p, n, b):
    """
    Computes Q given vf, vr, Efield, p, n, b, ro, and g.
    ***Chat GPT translated this code from mathematica to python for me. I'm not sure I trust it yet so be suspicious of it. -Lucas
    
    Parameters:
        vf (tuple): (value, uncertainty) for vf
        vr (tuple): (value, uncertainty) for vr
        Efield (tuple): (value, uncertainty) for Efield
        p (tuple): (value, uncertainty) for p
        n (tuple): (value, uncertainty) for n
        b (float): constant
    
    Returns:
        tuple: (propagated value, propagated uncertainty)
    """
    # Extract values and uncertainties
    vf_val, sigma_vf = vf
    vr_val, sigma_vr = vr
    Efield_val, sigma_Efield = Efield
    p_val, sigma_p = p
    n_val, sigma_n = n
    ro_val = 886
    g_val = 9.8

    # Compute A(n)
    A_n = np.sqrt((9 * n_val**3) / (2 * ro_val * g_val))
    # Compute function value Q
    prefactor = (6 * np.pi / Efield_val) * (vf_val + vr_val) * np.sqrt(vf_val)
    denominator = 1 + (b / (p_val * A_n))
    Q_val = prefactor * denominator**(-3/2) * A_n

    # Compute uncertainty propagation (translated from Mathematica)
    sigma_Q = np.sqrt(
        ((9 * np.sqrt(2) * np.pi * np.sqrt(vf_val) * np.sqrt(n_val**3 / (g_val * ro_val)) 
          * (denominator**(-3/2))) / Efield_val 
         + (9 * np.pi * (vf_val + vr_val) * np.sqrt(n_val**3 / (g_val * ro_val)) 
            * (denominator**(-3/2))) / (np.sqrt(2) * np.sqrt(vf_val) * Efield_val))**2 * sigma_vf**2
        + (162 * np.pi**2 * vf_val * n_val**3 * sigma_vr**2) / (g_val * ro_val * Efield_val**2 * denominator**3)
        + (162 * np.pi**2 * vf_val * (vf_val + vr_val)**2 * n_val**3 * sigma_Efield**2) 
          / (g_val * ro_val * Efield_val**4 * denominator**3)
        + (81 * b**2 * np.pi**2 * vf_val * (vf_val + vr_val)**2 * sigma_p**2) 
          / (Efield_val**2 * p_val**4 * denominator**5)
        + ((27 * np.pi * np.sqrt(vf_val) * (vf_val + vr_val) * n_val**2 * sigma_n**2 * denominator**(-3/2)) 
           / (np.sqrt(2) * g_val * ro_val * Efield_val * np.sqrt(n_val**3 / (g_val * ro_val))) 
           + (27 * b * np.pi * np.sqrt(vf_val) * (vf_val + vr_val) * denominator**(-5/2)) 
           / (2 * Efield_val * p_val * n_val))**2 * sigma_n**2
    )

    return Q_val, sigma_Q


def viscosity(T):
    """
    Computes the viscosity and its propagated uncertainty.
    
    Parameters:
    T (tuple): (value of temperature, uncertainty in temperature)
    C (float): Constant parameter C
    
    Returns:
    tuple: (viscosity value, propagated uncertainty)
    """
    T_val, sigma_T = T
    
    # Compute viscosity value
    viscosity_val = (1.827e-5) * ((291.15 + 120) / (T_val + 120)) * (T_val / 291.15) ** (3/2)
    
    # Compute partial derivative d(viscosity)/dT
    term1 = (-1.51204e-6 * T_val ** (3/2)) / (120 + T_val) ** 2
    term2 = (2.26806e-6 * np.sqrt(T_val)) / (120 + T_val)
    d_viscosity_dT = term1 + term2
    
    # Propagate uncertainty
    sigma_viscosity = np.sqrt((d_viscosity_dT * sigma_T) ** 2)
    
    return viscosity_val, sigma_viscosity

import numpy as np

def efield(V, d):
    """
    Computes the electric field and its propagated uncertainty.
    
    Parameters:
    V (tuple): (value of voltage, uncertainty in voltage)
    d (tuple): (value of distance, uncertainty in distance)
    
    Returns:
    tuple: (electric field value, propagated uncertainty)
    """
    V_val, sigma_V = V
    d_val, sigma_d = d
    
    # Compute electric field value
    E_val = V_val / d_val
    
    # Compute uncertainty propagation
    sigma_E = np.sqrt((sigma_V / d_val) ** 2 + ((V_val * sigma_d) / d_val ** 2) ** 2)
    
    return E_val, sigma_E

def resistance_to_temperature(resistance, resistance_values=[3.239,3.118,3.004,2.897,2.795,2.700,2.610,2.526,2.446,2.371,2.300,2.233,2.169],
                            temperature_values = [10,11,12,13,14,15,16,17,18,19,20,21,22]):
    '''
    # Example usage:
    Temp = [10,11,12,13,14,15,16,17,18,19,20,21,22]
    Resistance = [3.239,3.118,3.004,2.897,2.795,2.700,2.610,2.526,2.446,2.371,2.300,2.233,2.169]

    for i in range(len(r)):
        resistance_input = r[i]
        temperature_output = resistance_to_temperature(r[i], Resistance, Temp)
        print(f"Estimated temperature for BBL {str(i+1)}: {temperature_output:.2f}°C")
    
    BBL8 = BBL9
    BBL9 = BBL10
    BBL10 = BBL11
    BBL11 = BBL12
    '''
    temperature_values = np.array(temperature_values)+273.15
    return np.interp(resistance, resistance_values[::-1], temperature_values[::-1])


def inches_to_pascals(inches):
    # Conversion factor
    conversion_factor = 3386.39
    # Convert inches of mercury to pascals
    pascals = inches * conversion_factor
    return pascals

def velocity(y, t):
    """
    Computes the velocity and its propagated uncertainty.
    
    Parameters:
    y (tuple): (value of distance, uncertainty in distance)
    t (tuple): (value of time, uncertainty in time)
    *NOTE: uncertainty in time is assumed to be zero
    
    Returns:
    tuple: (velocity value, propagated uncertainty)
    """
    y_val, sig_y = y
    t_val, sig_t = t, 0
    
    # Compute velocity value
    v_val = y_val / t_val
    
    # Compute uncertainty propagation
    sigma_v = sig_y / t_val
    
    return v_val, sigma_v


def drop_charge_and_uncertainty(Efield, Vfall, Vrise, eta, pressure, rho=886, g=9.81, b=8.2*(10**-3) ):
    """
    Compute the drop-charge value Q and its uncertainty errQ.

    According to your latest formula, we consider uncertainties only in:
      - Vfall
      - Vrise
      - Efield
    and treat the uncertainties in b, eta, rho, g, pressure as negligible.

    Parameters
    ----------
    Efield : float or (float, float)
        Electric field, possibly (value, sigma_E)
    Vfall  : float or (float, float)
        Fall velocity, possibly (value, sigma_Vfall)
    Vrise  : float or (float, float)
        Rise velocity, possibly (value, sigma_Vrise)
    b      : float
        Correction constant b (assumed negligible uncertainty here)
    eta    : float
        Viscosity of the fluid (assumed negligible uncertainty)
    rho    : float
        Density of the fluid (assumed negligible uncertainty)
    g      : float
        Acceleration due to gravity (assumed negligible uncertainty)
    pressure : float
        Ambient pressure (assumed negligible uncertainty)

    Returns
    -------
    (Q, errQ) : (float, float)
        Q     = central drop-charge value
        errQ  = uncertainty in Q, based on partial derivatives w.r.t.
                Vfall, Vrise, and Efield only, per the provided formula.
    """

    # Helper to parse input as (value, sigma). If only a float is given,
    # we treat sigma=0.0.
    def parse_input(x):
        if isinstance(x, tuple):
            return x[0], x[1]
        else:
            return x, 0.0

    # Parse only the three that have non-negligible uncertainties:
    E_val,   E_sig   = parse_input(Efield)
    Vf_val,  Vf_sig  = parse_input(Vfall)
    Vr_val,  Vr_sig  = parse_input(Vrise)

    # The rest are floats with negligible uncertainty:
    b_val       = b
    eta_val     = eta
    rho_val     = rho
    g_val       = g
    pressure_val= pressure

    #--------------------------------------------------------------------------
    # 1) Compute central value of Q.
    #
    #    Q = (6 π / E) * (Vfall + Vrise) * sqrt(Vfall)
    #         * sqrt( 9 η^3 / (2 ρ g) )
    #         * [ 1 / ( 1 + b / (P sqrt(9 Vfall η / (2 ρ g))) ) ]^(3/2)
    #--------------------------------------------------------------------------
    numerator  = 6.0 * math.pi / E_val
    factor1    = (Vf_val + Vr_val)
    factor2    = math.sqrt(Vf_val)
    factor3    = math.sqrt(9.0 * (eta_val**3) / (2.0 * rho_val * g_val))
    denom_inner= 1.0 + b_val / (
        pressure_val * math.sqrt(9.0 * Vf_val * eta_val / (2.0 * rho_val * g_val))
    )
    factor4    = (1.0 / denom_inner)**1.5

    Q = numerator * factor1 * factor2 * factor3 * factor4

    #--------------------------------------------------------------------------
    # 2) Compute uncertainty errQ based on your new (shorter) formula:
    #
    #   errQ = sqrt(
    #      [ (∂Q/∂Vfall) ]^2 * (errVfall)^2
    #    + [ (∂Q/∂Vrise) ]^2 * (errVrise)^2
    #    + [ (∂Q/∂Efield)]^2 * (errEfield)^2
    #   )
    #
    # where the partial derivatives are exactly as in your provided expression:
    #
    #   errQ = sqrt(  [ ( ... )^2 * errVfall^2 ]
    #               + [ ( ... )^2 * errVrise^2 ]
    #               + [ ( ... )^2 * errEfield^2 ]  )
    #
    # and the big bracketed expressions are the partial derivatives from
    # your Mathematica snippet. We'll just replicate them carefully.
    #--------------------------------------------------------------------------
    # Make short references for readability:
    Vf   = Vf_val
    Vr   = Vr_val
    E    = E_val
    b_   = b_val
    et   = eta_val
    rh   = rho_val
    gg   = g_val
    P    = pressure_val

    # Common sub-expressions to reduce repetition:
    sqrtVf = math.sqrt(Vf)
    big_sqrt_factor = math.sqrt((et**3) * (1.0/(
        1.0 + (math.sqrt(2)*b_)/(3.0 * P * math.sqrt((Vf * et)/(gg * rh)))
    ))**1.5 / (gg * rh))

    # partial wrt Vfall (the big bracket in front of errVfall^2):
    partial_vfall = (
        ( 9*math.sqrt(2)*math.pi * sqrtVf * big_sqrt_factor )/E
        + (
            9*math.pi*(Vf + Vr)*big_sqrt_factor
        )/( math.sqrt(2)*E*math.sqrt(Vf) )
        + (
            9*b_*math.pi * math.sqrt(Vf)*(Vf + Vr)* (et**4)
            *(1.0/(1.0 + (math.sqrt(2)*b_)/(
                3.0*P*math.sqrt((Vf*et)/(gg*rh)))
             ))**2.5
        ) / (
            4.0*E*(gg**2)*P*((Vf*et)/(gg*rh))**1.5 *
            math.sqrt((et**3)*(1.0/(
                1.0 + (math.sqrt(2)*b_)/(
                    3.0*P*math.sqrt((Vf*et)/(gg*rh))
                )
            ))**1.5 / (gg*rh)) * (rh**2)
        )
    )

    # partial wrt Vrise (the next bracket):
    partial_vrise = (
        ( 9*math.sqrt(2)*math.pi * math.sqrt(Vf) * big_sqrt_factor )/ E
    )

    # partial wrt Efield (the last bracket):
    partial_efield = (
       -(
         ( 9*math.sqrt(2)*math.pi * math.sqrt(Vf) * (Vf + Vr) * big_sqrt_factor )
         /( E**2 )
       )
    )

    # Now assemble them into the error formula:
    errQ = math.sqrt(
        (partial_vfall  ** 2) * (Vf_sig**2) +
        (partial_vrise  ** 2) * (Vr_sig**2) +
        (partial_efield ** 2) * (E_sig **2)
        # If you later want to account for b_sig, just add a + (...)^2 * b_sig^2
    )

    return Q, errQ
