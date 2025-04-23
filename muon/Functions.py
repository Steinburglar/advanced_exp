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