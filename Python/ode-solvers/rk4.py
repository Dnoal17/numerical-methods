import numpy as np
#=================================================================
# RUNGE-KUTTA 4 ALGORITHM
#=================================================================
# Author: Daniel Noal Pineda
# Email : noaldaniel41@gmail.com
# Date  : 2026
# Repository: https://github.com/Dnoal17/numerical-methods.git
#=================================================================
# OBJECTIVES: Implement the 4th order Runge-Kutta method (RK4) 
#             to solve differential equations
#
# ERROR: O(dt^4) locally, O(dt^3) globally
#
# REQUIREMENTS : External function containing the system of
#                diferential equations to be solved is required
#
# INPUTS:
#         ·t    : current time
#         ·dt   : time step
#         ·y_in : input state as a numpy array 
#
# OUTPUTS:
#         ·y_out: output state as a numpy array after time dt
#=================================================================

import numpy as np

def rk4_step(t, dt, y_in, derivad):

    # Calculate the K_i of the RK4 algorithm using an auxiliary vector
    # Requires an external subroutine containing the system of differential equations

    # K1
    K1 = derivad(t, y_in)

    # K2
    yaux = y_in + 0.5 * dt * K1
    K2 = derivad(t + 0.5 * dt, yaux)

    # K3
    yaux = y_in + 0.5 * dt * K2
    K3 = derivad(t + 0.5 * dt, yaux)

    # K4
    yaux = y_in + dt * K3
    K4 = derivad(t + dt, yaux)

    # RK4 ALGORITHM
    y_out = y_in + (dt / 6.0) * (K1 + 2 * K2 + 2 * K3 + K4)

    return y_out