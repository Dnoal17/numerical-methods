#=================================================================
# 0 ORDER INTERPOLATION
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2026
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Interpolating f(x) at order 0 from a set of points 
#             {(xi,fi)}
#                       
# INPUTS:
#         ·X,Fi : two arrays containing the points such that f(Xi) = Fi, where Xi must be sorted
#         ·x    : the point where f(x) is to be interpolated
#
# OUTPUTS: 
#         ·f(x)   : 0 order interpolation of f at point x
#=================================================================

import numpy as np

def zero_order_interpolation(Xi, Fi, x):

    # Make sure Xi and Fi numpy arrays and inicialiez fx
    Xi = np.array(Xi)
    Fi = np.array(Fi)

    # Check both arrays have the same dimension
    if len(Xi) == len(Fi):

        npts = len(Xi)
        # Error message if x is not within the interval of points being considered
        if x < Xi[0] or Xi[npts - 1] < x:

            print('x is not within the interval considered')

            return

        else:

            # Find withing which two points does x lie
            for i in range(0, npts - 1):  

                if Xi[i] <= x <= Xi[i+1]:

                    return Fi[i]
                
    # Error message when dimensions are not equal
    else:

        print('The two arrays do not match in size')

        return
