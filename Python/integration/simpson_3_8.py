#=================================================================
# SIMPSON 3/8 ALGORITHM
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Integrating f(x) in the interval [a,b] using the 
#             Simpson 3/8 algorithm
#
# ERROR: O(h^4), with h the length of the subintervals 
#                     
# INPUTS:
#         ·a,b      : interval endpoints
#         ·k        : integer used to calculate the number of subintervals used (N=3*2^k)
#         ·fun      : 1 variable - callable function (NumPy-compatible)
#
# OUTPUTS: 
#         ·integral : integral aproximation
#=================================================================

import numpy as np

def simpson_3_8(a, b, k, fun):

    # Initialize the variables
    N = 3 * 2**k
    h = (b - a) / N
    X = np.linspace(a, b, N + 1)
    Y = fun(X)

    # Integral approximation
    idx = np.arange(1,N)
    integral = (3*h / 8) * (Y[0] + Y[-1] + 3 * np.sum(Y[idx % 3 != 0]) + 2 * np.sum(Y[idx % 3 == 0]))

    return integral