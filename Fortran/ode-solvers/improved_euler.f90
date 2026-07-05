!=================================================================
! IMPROVED EULER METHOD ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Apply Improved Euler's method to solve differential equations 
!             as a system of equations
!
! ERROR: O(h), with h the step size
!
! INPUTS:
!         ·N    : number of points used
!         ·y    : initial conditions array (N x nvar)
!         ·fun  : external Fortran function containing the system of diferential equations to be solved
!         ·a,b  : interval endpoints
!         ·nvar : number of variables in the system
!
! OUTPUTS:
!         ·y    : solution array (updated with the numerical solution)
!=================================================================

subroutine imporved_euler_method(N, y, fun, a, b, nvar)

    implicit none

    ! Inputs
    integer, intent(in) :: N, nvar
    double precision, intent(in) :: a, b
    double precision, intent(inout) :: y(N, nvar)
    external :: fun

    ! Internal Variables
    double precision :: h, t
    integer :: i, j

    ! Determine the step h
    h = (b - a)/dble(N - 1)

    ! Use the basic Euler algorithm to initialize y(2,:), necessary for the improved algorithm
    do j = 1, nvar
        y(2,j) = y(1,j) + h * fun(a, y(1,:), j)
    end do

    ! Apply the improved Euler method algorithm. Start the loop at i=2 because we start from initial conditions
    do i = 2, N-1

        ! Calculate the value of the independent variable for each iteration
        t = a + h*dble(i-1)

        ! Update each variable of the differential equation system
        ! The variable j serves to change the equation within the system, allowing second order equations to be solved (which are the ones we are interested in)
        do j = 1, nvar

            y(i+1,j) = y(i-1,j) + 2.0d0 * h * fun(t, y(i,:), j)

        end do

    end do

    return
    
end subroutine imporved_euler_method