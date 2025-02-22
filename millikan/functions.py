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
        return (np.mean(data), np.std(data, ddof=1))  # ddof=1 for sample std dev

    # Compute mean and standard deviation for rising and falling droplets
    rise_stats = compute_mean_stddev(v_rise_trimmed)
    fall_stats = compute_mean_stddev(v_fall_trimmed)

    return rise_stats, fall_stats