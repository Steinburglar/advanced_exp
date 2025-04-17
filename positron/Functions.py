"""contains functions for the positron analysis
"""
import numpy as np
import math


def circles_overlap(angle, r1, r2, length):
    """returns the overlapping area of two circles of radius r1 and r2, wth a distance d between their centers
    """
    d = d_from_angle(angle, length)
    if d>= r1 + r2:
        #print(f"No overlap at {angle} degrees")
        return 0
    elif d <= np.abs(r1 - r2):
        #print(f"One circle is inside the other at {angle} degrees")
        return math.pi * min(r1, r2)**2
    else:
        d1 = (r1**2 - r2**2 + d**2)/(2*d)
        d2 = d - d1
        term1 = r1**2 * math.acos(d1/r1)
        term2 = r2**2 * math.acos(d2/r2)
        term3 = d1 * math.sqrt(r1**2 - d1**2)
        term4 = d2 * math.sqrt(r2**2 - d2**2)
        area = term1 + term2 - term3 - term4
        #print(f"Overlap at {angle} degrees: {area}")
        return area

def d_from_angle(angle, length):
    """returns the distance between centers of apperature circles when projected on a flat plane

    Args:
        angle (_float_): angle in degrees of arm
        length (_type_): true length from pivot to apperature, i.e not to opening of columnator but to the detector itself
    """
    radians = np.abs(math.radians(angle))
    d = length * math.tan(radians)
    return d

def normalized_overlap(angle, r1=1.75, r2=0.6, length=50):
    """returns the normalized overlapping area of two circles of radius r1 and r2, wth a distance d between their centers
    """
    area = circles_overlap(angle, r1, r2, length)
    norm_area = area / (math.pi * min(r1, r2)**2)
    return norm_area

