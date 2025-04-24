import numpy as np

def exponential_decay(t, a, tau):
    """
    Exponential decay function
    :param t: time value
    :param a: initial amplitude
    :param tau: decay constant
    :return: y value
    """
    return a * np.exp(-t/tau)

def corrected_exponential_decay(t, a, tau, gamma):
    """
    Corrected exponential decay function
    :param t: time value
    :param a: initial amplitude
    :param tau: decay constant
    :return: y value
    """
    return a * np.exp(-t/tau) + gamma*t