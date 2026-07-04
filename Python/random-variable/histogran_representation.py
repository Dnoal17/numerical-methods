#=================================================================
# HISTOGRAM REPRESENTATION
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Constructing a histogram from a data set and estimating 
#             the statistical uncertainty of each bin
#
# INPUTS:
#         ·x_data    : numpy array containing the data sample
#         ·n_box     : number of histogram bins
#
# OUTPUTS:
#         ·xhisto   : center of each histogram bin
#         ·histo    : normalized histogram (probability density)
#         ·errhisto : estimated error of each bin
#=================================================================

import numpy as np

def histogram(x_data, n_box):

    # Determine the minimum and maximum values of the data sample
    x_min = np.min(x_data)
    x_max = np.max(x_data)

    # Define the bin width using the selected number of bins
    h = (x_max - x_min) / n_box

    # Initialize all output arrays
    x_histo = np.linspace(x_min + h/2, x_max - h/2, n_box)    
    histo = np.zeros(n_box)
    err_histo = np.zeros(n_box)

    # Calculate the number of points for further use
    n_data = x_data.size

    # Identify the bin corresponding to each data point
    for i in range(n_data):

        # The maximum value must be treated separately to avoid assigning it to a non-existent bin
        if x_data[i] == x_max:
            i_box = n_box - 1

        # General case
        else:
            i_box = int((x_data[i] - x_min) / h)

        # Increase the counter of the corresponding bin
        histo[i_box] += 1

    # Normalize the histogram to obtain a probability density
    histo = histo / (n_data * h)

    # Estimate the statistical uncertainty assuming binomial statistics in each bin
    p_i = histo * h
    err_histo = np.sqrt(p_i * (1.0 - p_i) / n_data) / h

    return x_histo, histo, err_histo