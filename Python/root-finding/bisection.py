#=================================================================
# BISECTION METHOD
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Finding a root of f(x) = 0 in the interval [A,B] 
#             using the bisection method
#
# REQUIREMENTS: f(A)·f(B) < 0
#
# CONVERGENCE ORDER: The error in this method is bounded by the
#                    number of iterations as: error <= (B-A)/2^n
#                       
# INPUTS:
#         ·A,B      : the points that define the integration interval
#         ·fun      : defined as an external function [fun(x,f(x)]
#         ·er_bound : desired error bound
#
# OUTPUTS: 
#         ·root     : estimated root within the error bound
#         ·n_iter    : number of iterations needed
#=================================================================
import numpy as np

def bisection(A, B, er_bound, fun):

    # Inicialize the values of the images at the endpoints of the interval
    fA = fun(A)
    fB = fun(B)

    # Calculate the number of iterations needed to achieve the error bound
    n_iter = int(np.log((B - A) / er_bound) / np.log(2)) + 1

    # Check the sign requirement is fulfilled
    if (fA * fB < 0):
        
        for i in range (n_iter + 1):
            
            # Calculate the midpoint, C, and the image
            C = (A + B) / 2
            fC = fun(C)

            # Find which new interval, [A,C] or [C,B], is the root in by checking where there is a change of sign

            if (fC * fA < 0.0):

                # Prepare the next iterations
                B = C
                fB = fC

            elif (fC * fB < 0.0):

                # Prepare the next iterations
                A = C
                fA = fC

        # Take the last midpoint as the root estimation
        return C
    
    # If the input does not fulfill the requirement notify the reasson why 
    elif (fA * fB > 0):

        print(f"ERROR: Both A = {A} and B = {B} have the same sign.")

    # It may be the case that A or B are the root themselves
    elif (fA == 0.0):

        return A

    elif (fB ==0.0):
        
        return B