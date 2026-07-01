
#=================================================================
# TRAPEZOIDAL ALGORITHM
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2026
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Integrating f(x) in the interval [a,b] using the
#             trapezoidal rule
#
# ERROR: O(h^2), where h is the step size
#
# INPUTS:
#         ·a, b   : interval endpoints
#         ·k      : determines number of subintervals (N = 2^k)
#         ·fun    : 1 variable - callable function (vectorized, NumPy-compatible)
#
# OUTPUT:
#         ·integral : approximation of the integral
#=================================================================

import numpy as np


def trapezoidal(a, b, k, fun):

    N = 2**k
    h = (b - a) / N

    x = np.linspace(a, b, N + 1)
    y = fun(x)

    integral = h * (0.5 * (y[0] + y[-1]) + np.sum(y[1:-1]))

    return integral