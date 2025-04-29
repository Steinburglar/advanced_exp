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

def area_overlap(angle, r1=2.35, r2=0.95, length_1=44.45, length_2=23.9):
    """returns the non-normalized overlapping area of two circles of radius r1 and r2, wth a distance d between their centers
    """
    ratio = length_1 / length_2
    r2_scaled = 1.57#r2 * ratio
    area = circles_overlap(angle, r1, r2_scaled, length_1,)
    return area

