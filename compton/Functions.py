import pandas as pd
import math
import numpy as np
from cavendish.utils.Functions import *
from millikan.functions import *
from scipy.special import erf

def gaussian(x, a, mu, sigma):
    return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def modified_gaussian(x, a, mu, sigma, b):
    return (a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2)) )+ b

def compton_edge(E_gamma, m_e=511):
    return E_gamma / (1 + (E_gamma / (2 * m_e)))

def safe_divide(numerator, denominator):
    """
    Safely divides the numerator by the denominator.
    If the denominator is 0, it returns the numerator divided by 1 instead.

    Args:
        numerator (float): The numerator of the division.
        denominator (float): The denominator of the division.

    Returns:
        float: The result of the division.
    """
    if denominator == 0:
        return numerator / 1
    else:
        return numerator / denominator
    
def logistic(x, a, b, c, d):
    return a / (1 + np.exp(-c * (x - d))) + b

def linear(x, m, b):
    return m * x + b

def fit_erf(x, a, b, c, d):
    return a * erf((x - d) / c) + b