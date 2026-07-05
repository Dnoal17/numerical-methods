!=================================================================
! TRAPEZOIDAL ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Integrating f(x) in the interval [a,b] using the
!             Trapezoidal algorithm
!
! ERROR: O(h^2), with h the length of the subintervals
!
! INPUTS:
!         ·a,b      : interval endpoints
!         ·k        : integer used to calculate the number of subintervals used (N=2^k)
!         ·fun      : defined as an external Fortran function [fun(x)]
!
! OUTPUTS: 
!         ·integral : integral aproximation
!=================================================================

subroutine trapezoidal(a, b, k, fun, integral)

    implicit none

    ! Inputs
    double precision, intent(in) :: a,b
    integer, intent(in) :: k
    double precision, external :: fun

    ! Outputs
    double precision, intent(out) :: integral

    ! Internal Variables
    double precision :: h, x, sum
    integer :: i, N

    ! Initialize the variables
    sum = 0.0d0
    N = 2**k
    h = (b-a)/dble(N)

    ! Iterate through the number of subintervals
    do i=0,N

        ! Define the points
        x = a + h*dble(i)

        ! TRAPEZOIDAL ALGORITHM

        ! Contribution of the endpoints
        if (i == 0 .or. i == N) then
            sum = sum + 0.5d0*fun(x)

        ! Contribution of interior points
        else
            sum = sum + fun(x)
        endif

    enddo

    ! Integral estimation
    integral = sum*h

    return

end subroutine TRAPEZOIDAL