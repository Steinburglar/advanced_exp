"""This script contains functions for use in other parts of the analysis.
"""

import numpy as np

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
