# Description: Auxiliary functions for script 'precipitation_indices.py'
#
# Author: Christian R. Steger, March 2023

# Load modules
import numpy as np
import numba as nb


# -----------------------------------------------------------------------------

@nb.jit((nb.float32[:, :, :], nb.float32[:, :, :], nb.int64, nb.int64,
         nb.int64), nopython=True, parallel=True)
def update_max_values_all_day(prec_keep, prec, len_y, len_x, num_keep):
    """Update maximal precipitation values for all day/hour percentile
    calculation.

    Parameters
    ----------
    prec_keep : ndarray of float
        Array (three-dimensional) with retained precipitation data (y, x, time)
    prec : ndarray of float
        Array (three-dimensional) with precipitation data (time, y, x)
    len_y : int
        Dimension length in y-direction
    len_x : int
        Dimension length in x-direction
    num_keep : int
        Number of elements to keep"""

    for i in nb.prange(len_y):
        for j in range(len_x):
            mask = (prec[:, i, j] > prec_keep[i, j, 0])
            prec_keep[i, j, :] \
                = np.sort(np.append(prec_keep[i, j, :],
                                    prec[mask, i, j]))[-num_keep:]


# -----------------------------------------------------------------------------

@nb.jit((nb.float32[:, :, :], nb.float32[:, :, :], nb.int64, nb.int64,
         nb.int64, nb.float64), nopython=True, parallel=True)
def update_max_values_wet_day(prec_keep, prec, len_y, len_x, num_keep,
                              prec_thresh):
    """Update maximal precipitation values for wet day/hour percentile
    calculation.

    Parameters
    ----------
    prec_keep : ndarray of float
        Array (three-dimensional) with retained precipitation data (y, x, time)
    prec : ndarray of float
        Array (three-dimensional) with precipitation data (time, y, x)
    len_y : int
        Dimension length in y-direction
    len_x : int
        Dimension length in x-direction
    num_keep : int
        Number of elements to keep
    prec_thresh : float
        Threshold for wet day/hour"""

    for i in nb.prange(len_y):
        for j in range(len_x):
            mask = (prec[:, i, j] > prec_keep[i, j, 0]) \
                   & (prec[:, i, j] > prec_thresh)
            prec_keep[i, j, :] \
                = np.sort(np.append(prec_keep[i, j, :],
                                    prec[mask, i, j]))[-num_keep:]


# -----------------------------------------------------------------------------

@nb.jit(nopython=True)
def weighted_average(data_in, weights, weight_zero="set_nan"):
    """Compute weighted average over first axis. Weights of 0 are treated
    according to the argument 'weight_zero'.

    Parameters
    ----------
    data_in : ndarray of float
        Array (three-dimensional) with input data (time, y, x)
    weights : ndarray of int
        Array (three-dimensional) with weights (time, y, x)
    weight_zero : str
        Rule for weights of 0. For 'set_nan', any weights equal to 0 will yield
        NaN as an output. For 'ignore', weights equal to 0 are ignored. In case
        all weights are 0, the output will be NaN.

    Returns
    -------
    data_out : ndarray of float
        Array (two-dimensional) with weighted output data (y, x)"""

    # Check arguments
    if weight_zero not in ("set_nan", "ignore"):
        raise ValueError("Invalid argument for 'weight_zero'")

    # Compute weighted average
    data_out = np.empty(data_in.shape[1:], dtype=np.float32)
    data_out[:] = np.nan
    if weight_zero == "set_nan":
        for i in range(data_in.shape[1]):
            for j in range(data_in.shape[2]):
                if np.any(weights[:, i, j] == 0):
                    continue
                weights_norm = weights[:, i, j] / np.sum(weights[:, i, j])
                data_out[i, j] = np.sum(data_in[:, i, j] * weights_norm)
    else:
        for i in range(data_in.shape[1]):
            for j in range(data_in.shape[2]):
                if np.all(weights[:, i, j] == 0):
                    continue
                mask = weights[:, i, j] > 0
                weights_norm = weights[mask, i, j] \
                    / np.sum(weights[mask, i, j])
                data_out[i, j] = np.sum(data_in[mask, i, j] * weights_norm)

    return data_out
