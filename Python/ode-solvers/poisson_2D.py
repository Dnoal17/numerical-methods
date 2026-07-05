#=================================================================
# POISSON EQUATION SOLVER ALGORITHM
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2026
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Solve the Heat Equation (2D Poisson equation) using Jacobi, 
#             Gauss-Seidel, and Successive Over-Relaxation (SOR) methods
#
# NOTES: The method used is controlled by the variable icontrol
# 
# INPUTS:
#         ·icontrol     : method selector (1=Jacobi, 2=Gauss-Seidel, 3=SOR)
#         ·t_cinitials  : initial temperature matrix
#         ·p            : source term matrix
#         ·h            : grid spacing
#         ·i_point      : i-index for convergence monitoring point
#         ·j_point      : j-index for convergence monitoring point
#         ·omega        : relaxation parameter (for SOR method)
#         ·unit_file    : file unit for writing convergence data
#
# OUTPUTS:
#         ·T_final      : final temperature matrix (solution)
#         ·conv_p       : convergence vector at the monitoring point
#=================================================================

import numpy as np

def solver_poisson_2d(icontrol, t_cinitials, p, h, i_point, j_point, omega):

    # Convergence parameters
    precision = 1.0e-4
    max_iter = 100000

    # Matrix dimensions
    Nx, Ny = t_cinitials.shape

    # Initialize convergence vector
    conv_p = np.zeros(max_iter)

    # Initialization
    convergence = False
    iteration = 0

    # Copy initial conditions
    T_k = t_cinitials.copy()
    T_k_plus_1 = T_k.copy()

    # Iterate until convergence
    while not convergence and iteration < max_iter:

        iteration += 1

        #-------------------------------
        # JACOBI METHOD (ICONTROL = 1)
        #-------------------------------

        if icontrol == 1:

            # Jacobi method algorithm
            T_k_plus_1[1:-1, 1:-1] = (T_k[2:, 1:-1] + T_k[:-2, 1:-1] + T_k[1:-1, 2:] +
                T_k[1:-1, :-2] +  h**2 * p[1:-1, 1:-1]) / 4.0

        #-------------------------------------
        # GAUSS-SEIDEL METHOD (ICONTROL = 2)
        #-------------------------------------

        elif icontrol == 2:

            for i in range(1, Nx - 1):
                for j in range(1, Ny - 1):

                    # Gauss-Seidel algorithm
                    T_k_plus_1[i, j] = (T_k_plus_1[i - 1, j] + T_k[i + 1, j] +
                        T_k_plus_1[i, j - 1] + T_k[i, j + 1] + h**2 * p[i, j]) / 4.0

        #---------------------------------------
        # SOR METHOD (ICONTROL = 3)
        #---------------------------------------

        elif icontrol == 3:

            for i in range(1, Nx - 1):
                for j in range(1, Ny - 1):

                    # Successive Over-Relaxation (SOR) algorithm
                    T_k_plus_1[i, j] = T_k[i, j] + omega * (T_k_plus_1[i - 1, j] +
                        T_k[i + 1, j] + T_k_plus_1[i, j - 1] +T_k[i, j + 1] - 4.0 * T_k[i, j] +
                        h**2 * p[i, j]) / 4.0
        
        # Error if icontrol != 1,2,3
        else:

            print(f"ERROR: Icontrol = {icontrol}. Its value must be 1, 2, 3")
            return None

        #------------------------------------------------
        # CONVERGENCE - ERROR CALCULATION - NEXT ITERATION
        #------------------------------------------------

        # Save to convergence vector
        conv_p[iteration - 1] = T_k_plus_1[i_point, j_point]

        # Calculate the error
        error = np.max(np.abs(T_k_plus_1[1:-1, 1:-1] - T_k[1:-1, 1:-1]))

        # Check convergence
        if error < precision:
            convergence = True

        # Prepare the next iteration
        T_k = T_k_plus_1.copy()

    # Save the final result
    T_final = T_k_plus_1.copy()

    return T_final, conv_p[:iteration]