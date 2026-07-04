#=================================================================
# NEWTON-RAPHSON METHOD
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2025
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Finding a root of f(x) = 0 using the Newton-Raphson
#             algorithm
#
# REQUIREMENTS: f'(x0) not being too close to zero
#
# CONVERGENCE ORDER: Quadratic convergence near the root
#                       
# INPUTS:
#         ·x0       : the point where the iteration begins
#         ·fun      : 1D variable - callable function (NumPy-compatible) with output (f(x),f'(x)) in a tuple
#         ·max_iter : maximum number of iterations permited
#         ·er_bound : desired error bound
#
# OUTPUTS: 
#         ·root     : estimated root within the error range
#         ·n_iter    : number of iterations needed
#=================================================================

def newtonrap(x0, fun, max_iter, er_bound):

    # Initialize x with the given value
    x = x0
    
    for n_iter in range(max_iter): 

        # Calculate the values of f(x) and f'(x)
        fx, dfx = fun(x)
        
        # If f'(x) is too small the alorithm will diverge and we must stop the iterations
        if abs(dfx) < 1e-13:

            print("DIVERGENCE ERROR: f'(x) too small")

            return None
        
        # NEWTON-RAPHSON ALGORITHM
        root = x - fx / dfx
        
        # Stop the program once within the wished error bound
        if abs(root - x) < er_bound:

            return root, n_iter + 1
        
        # Prepare the next iteration
        x = root
    
    # Error code if the maximum number of iterations is reached
    print("ERROR: Maximum number of iterations reached")
    return None

