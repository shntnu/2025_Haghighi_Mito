"""Vectorized implementation of slope calculation for radial distribution analysis.

This module provides both the original row-by-row implementation and a vectorized
version for batch processing of radial distribution profiles. The vectorized version
provides significant speedup when processing large datasets.
"""

import numpy as np
from scipy.signal import savgol_filter


def smooth_data(data, window_length=5, polyorder=3):
    """Apply Savitzky-Golay filter to smooth data.

    Parameters
    ----------
    data : array_like
        Input data to smooth
    window_length : int, optional
        Length of the filter window (must be odd), by default 5
    polyorder : int, optional
        Order of the polynomial used to fit the samples, by default 3

    Returns
    -------
    ndarray
        Smoothed data
    """
    return savgol_filter(data, window_length, polyorder)


def find_end_slope2(data, height=None):
    """Calculate slope from last peak to end of radial distribution profile.

    This is the original implementation that processes one row at a time.

    Parameters
    ----------
    data : array_like
        1D array representing a radial distribution profile
    height : float, optional
        Peak height threshold (unused in current implementation)

    Returns
    -------
    tuple of (int, float)
        last_peak_ind : Index of the last peak in the profile
        slope : Calculated slope from last peak to average of last two points
    """
    data = smooth_data(data)
    min_max_indc = [np.argmax(data), np.argmin(data)]
    last_peak_ind0 = [i for i in min_max_indc if i < len(data) - 2]
    if last_peak_ind0 == []:
        return 0, 0
    last_peak_ind = np.max(last_peak_ind0)
    last_two_points_amplitude = (data[-1] + data[-2]) / 2
    slope = (last_two_points_amplitude - data[last_peak_ind]) / (len(data) - last_peak_ind - 1)

    return last_peak_ind, slope


def find_end_slope2_vectorized(data_matrix):
    """Vectorized version of find_end_slope2 for batch processing.

    Processes all rows at once for significant speedup compared to row-by-row
    processing with np.apply_along_axis.

    Parameters
    ----------
    data_matrix : ndarray
        2D array where each row is a radial distribution profile to process

    Returns
    -------
    ndarray
        2D array with shape (n_rows, 2) containing [last_peak_ind, slope] for each row

    Examples
    --------
    >>> data = np.random.rand(1000, 16)  # 1000 profiles with 16 bins
    >>> results = find_end_slope2_vectorized(data)
    >>> last_peak_indices = results[:, 0]
    >>> slopes = results[:, 1]
    """
    # Smooth all rows at once - VECTORIZED
    smoothed = savgol_filter(data_matrix, window_length=5, polyorder=3, axis=1)

    n_rows, n_cols = smoothed.shape

    # Find argmax and argmin for each row - VECTORIZED
    argmax_inds = np.argmax(smoothed, axis=1)
    argmin_inds = np.argmin(smoothed, axis=1)

    # Determine which peaks are valid (not in the last 2 positions)
    argmax_valid = argmax_inds < (n_cols - 2)
    argmin_valid = argmin_inds < (n_cols - 2)

    # For rows where both are valid, take the max
    # For rows where only one is valid, take that one
    # For rows where neither is valid, use 0
    last_peak_ind = np.where(
        argmax_valid & argmin_valid,
        np.maximum(argmax_inds, argmin_inds),
        np.where(argmax_valid, argmax_inds, np.where(argmin_valid, argmin_inds, 0)),
    )

    # Check if there are any valid peaks
    has_valid_peak = argmax_valid | argmin_valid

    # Calculate slopes - VECTORIZED
    last_two_points_amplitude = (smoothed[:, -1] + smoothed[:, -2]) / 2

    # Get the value at the last peak for each row using advanced indexing
    peak_values = smoothed[np.arange(n_rows), last_peak_ind]

    # Calculate slope for all rows at once
    slopes = (last_two_points_amplitude - peak_values) / (n_cols - last_peak_ind - 1)

    # Set slope to 0 where there's no valid peak
    slopes = np.where(has_valid_peak, slopes, 0)
    last_peak_ind = np.where(has_valid_peak, last_peak_ind, 0)

    # Stack results into [last_peak_ind, slope] pairs
    results = np.column_stack([last_peak_ind, slopes])

    return results
