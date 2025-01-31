"""
Lucas Steinberger
1/31/2025
Python script containing  functions to for loading and converting data from Casey, Lucas, and Cooke's Cavendish lab in PHYS 064.
Goal is to provide utility to load csv and return 4 arrays, each containing the data points with corresponding uncertainties for a measurement,
with position converted to radians of beam rotation. Note that that means the values returned are *half* of the angle produced by the laser, since the law of reflection
will double the angle between the laser rays. 
"""

def load_convert(path, sig_x = 0, sig_s = 0, arrays = 4 ):
    """ Main function for dataloader.

    Args:
        path (_str_): Path to the location of the .csv file containing measurements of time and position (x) in cm. 
        

    Returns:
        _array_: Array of data arrays:  [data_aray1, data_aray2, data_array3 ....]
            
            _data_array_:
        
    """
    
    
    a = 2
    return a