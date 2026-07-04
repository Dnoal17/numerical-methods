#=================================================================
# BOX-MULLER METHOD
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Generate Gaussian random numbers N(0,1) using the 
#             Box-Müller transformation method
#
# NOTES: To generate Gaussian numbers with other constants sigma, mu:
#        N(mu,sigma) = sigma*N(0,1) + mu
#
# INPUTS:
#         ·ndat  : number of points to generate
#
# OUTPUTS:
#         ·xgaus : Numoy array containing the generated Gaussian numbers N(0,1)
#=================================================================

import numpy as np

def gauss(ndat):

    # The algorithm need ndat to be even. In case the user sets an odd number one extra point will be generated
    if ndat % 2 != 0:
        m = ndat + 1

    else:
        m = ndat

    # Generate m random numbers in U(0,1)
    u = np.random.rand(m)

    # Determine the parameters of the Box-Muller method, r and theta, and initialize the output array
    r = np.sqrt(-2 * np.log(u[0::2]))
    theta = 2 * np.pi * u[1::2]
    x = np.empty(m)

    # Generate the gaussian numbers N(0,1)
    x[0::2] = r * np.cos(theta)
    x[1::2] = r * np.sin(theta)

    # Return the array with the length the user set
    return x[:ndat]