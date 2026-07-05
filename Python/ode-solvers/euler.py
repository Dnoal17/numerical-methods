#=================================================================
# EULER METHOD
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Apply Euler's method to solve differential equations 
#             as a system of equations
#
# ERROR: O(h), with h the step size
#
# INPUTS:
#         ·y_in  : initial conditions as a numpy array (N x nvar)
#         ·fun   : external numpy-compatible function containing the system of diferential equations to be solved
#         ·a,b   : interval endpoints
#
# OUTPUTS:
#         ·y_out : solution array (updated with the numerical solution)
#=================================================================

import numpy as np

def euler(y_in, a, b, fun):

    # Determine the number of points used and step
    N = y_in.shape[0]
    h = (b - a) / (N - 1)

    # Inicialize the independent variable array
    t = np.linspace(a, b, N)

    # Initialize output array
    y_out = y_in.copy()

    # EULER ALGORITHM
    for i in range(1, N):
        y_out[i] = y_out[i-1] + h * fun(t[i-1], y_out[i-1])

    return y_out