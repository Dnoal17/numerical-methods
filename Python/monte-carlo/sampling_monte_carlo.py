#=================================================================
# MONTECARLO SAMPLING ALGORITHM
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Integrating f(x) in N-dimensional interval through 
#             the Monte Carlo algorithm with importance sampling
#             using a probability density function
#
# REQUIREMENT: Density fucntion must be normalized
#
# ERROR: The more similar the function and the probability density 
#        are, the more efficient the method will be.
#
# INPUTS:
#         ·fun     : n-D function to be integrated - callable function (NumPy-compatible)
#         ·density : n-D probability density function - callable function (NumPy-compatible)
#         ·nums    : matrix containing the random numbers (N * dim)
#
# OUTPUTS: 
#         ·integral : integral approximation
#         ·error    : error estimation
#=================================================================

import numpy as np

def sampling_montecarlo(nums, fun, density):
    
    # Evaluate the integrand and the probability density at the sampled points
    f = fun(nums)
    g = density(nums)
    w = f / g

    # Calculate the integral value sampling with g
    integral = np.mean(w)

    # Estimated error: sigma / sqrt(N)
    sigma = np.sqrt(np.mean(w**2) - np.mean(w)**2)
    error = sigma / np.sqrt(len(w))

    return integral, error