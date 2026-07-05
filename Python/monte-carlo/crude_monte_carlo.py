#=================================================================
# MONTECARLO ALGORITHM
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2026
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Integrating f(x) in N-Dimensions using the Monte Carlo 
#             algorithm
#
# ERROR: O(1/√N), with N the number of random points (N>>1)
#                     
# INPUTS:
#         ·a,b      : vnumpy arrays defining the endpoints of all n integrals
#         ·N        : number of random points used each integration(integeer)
#         ·fun      : n-D function to be integrated - callable function (NumPy-compatible)
#
# OUTPUTS: 
#         ·integral : integral aproximation
#         ·error    : error estimation
#=================================================================

import numpy as np

def crude_montecarlo(a, b, N, fun):

    # Generate N*dim random numbers distributed uniformly in each interval and find the images
    dim = a.size
    X = a + (b - a) * np.random.rand(N, dim)
    f = fun(X)
    
    # Calculate the volume of the domain 
    lambda_ = np.prod(b - a)

    # Integral estimation
    integral = lambda_ * np.mean(f)

    # Error estimation: lambda * sgima / √N
    sigma = np.sqrt(np.mean(f**2) - np.mean(f)**2)
    error = lambda_ * sigma / np.sqrt(N)

    return integral, error