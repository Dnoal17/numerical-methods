#=================================================================
# SECANT METHOD
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Finding a root of f(x) = 0 using the secant method.
#
# REQUIREMENTS: abs(f(x1)-f(x0)) being not to close to zero
#
# CONVERGENCE ORDER: Superlineal convergence near the root
#                       
# INPUTS:
#         ·x0,x1    : points where the iteration begins
#         ·fun      : 1D variable - callable function (NumPy-compatible) with output (f(x),f'(x)) in a tuple
#         ·max_iter : maximum number of iterations permited
#         ·er_bound : desired error bound
#
# OUTPUTS: 
#         ·root     : estimated root within the error range
#         ·niter    : number of iterations needed
#=================================================================


def secant(x0, x1, fun, max_iter, er_bound):

    # Initialize with the given points x0,x1
    x = x1
    x_prev = x0

    for n_iter in range(max_iter): 

        # Calculate the values of f(x) and f'(x)
        f1 = fun(x)
        f0 = fun(x_prev)

        # If f'(x) is too small the alorithm will diverge and we must stop the iterations
        if abs(f1-f0) < 1e-13:

            print("ERROR: denominator too small")
            return None
        
        # SECANT ALGORITHM
        root = x - f1 * (x - x_prev) / (f1 - f0)
        
        # Stop the program once within the wished error bound
        if abs(root - x) < er_bound:

            return root, n_iter + 1
        
        # Prepare the next iteration
        x_prev = x
        x = root

    # Error code if the maximum number of iterations is reached
    print("ERROR: Maximum number of iterations reached")
    return None