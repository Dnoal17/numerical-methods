#=================================================================
# MOMENT CALCULATOR
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2026
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Calculate the m-th order central moment, expected value 
#             and standard deviation of a set of random numbers 
#
# INPUTS:
#         ·x_data          : numpy array containing the data
#         ·m               : order of the central moment to calculate
#
# OUTPUTS:
#         ·expected_x      : expected value
#         ·moment_m        : m-th order central moment
#         ·stnd_deviation  : standard deviation
#=================================================================

import numpy as np

def moment_order_M(x_data, m):

    # 1) EXPECTED VALUE
    expected_x = np.mean(x_data)

    # 2) M-MOMENT
    moment_m = np.mean((x_data - expected_x)**m)

    # 3) Standard deviation
    stnd_deviation = np.sqrt(np.mean((x_data - expected_x)**2))

    return expected_x, moment_m, stnd_deviation