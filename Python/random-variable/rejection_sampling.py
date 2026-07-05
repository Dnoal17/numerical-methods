#=================================================================
# ACCEPTANCE-REJECTION METHOD
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2026
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Generate random numbers following a given probability 
#             density function g(x) using the acceptance-rejection method
#
# NOTES:
#         The probability density function must be defined (and normalized)
#         in the same interval [a,b] where the numbers will be generated
#
# INPUTS:
#         ·n_dat : number of points to generate
#         ·fun   : 1 variable - callable function f(x)
#         ·a,b   : interval endpoints
#         ·M     : upper bound such that fun(x) <= M for all x in [a,b]
#
# OUTPUTS:
#         ·xnums : array containing the generated random numbers
#=================================================================

import numpy as np

def subair(n_dat, fun, a, b, M):

    # Initialize the output array
    xnums = np.empty(n_dat)
    i = 0

    # Use a loop that doesn't stop until we have n_dat 
    while i < n_dat:

        # Generate 2 random numbers in U(0,1)
        x = np.random.rand()
        p = np.random.rand()

        # Convert these numbers to uniform distributions in [a,b] and [0,M]
        x = a + (b - a) * x
        p = M * p

        # Acceptance condition of the method
        if fun(x) >= p:

            xnums[i] = x

            # If we accept a number, count +1 in the loop
            i += 1

    return xnums