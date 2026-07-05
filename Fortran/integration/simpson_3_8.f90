!=================================================================
! SIMPSON 3/8 ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Integrating f(x) in the interval [a,b] using the 
!             Simpson 3/8 algorithm
!
! ERROR: O(h^4), with h the length of the subintervals 
!                     
! INPUTS:
!         ·a,b      : interval endpoints
!         ·k        : integer used to calculate the number of subintervals used (N=3*2^k)
!         ·f(x)     : defined as an external Fortran function [f(x)]
!
! OUTPUTS: 
!         ·integral : integral aproximation
!=================================================================

subroutine simpson_3_8(a, b, k, fun, integral)

    implicit none

    ! Inputs
    double precision, intent(in) :: a,b
    integer, intent(in) :: k
    external :: fun

    ! Outputs
    double precision, intent(out) :: integral

    ! Internal Variables
    double precision :: h, x, sum
    integer :: i, N

    ! Initialize the variables
    sum = 0.0d0
    N = 3*2**k
    h = (b-a)/dble(N)

    ! Iterate through the number of subintervals
    do i=0,N

        ! Define the points
        x = a + h*dble(i)

        ! SIMPSON 3/8 ALGORITHM

        ! Contribution of the endpoints
        if (i == 0 .or. i == N) then
            sum = sum + fun(x)

        ! Contribution of points labled by multiple of 3
        else if (mod(i,3) == 0) then
            sum = sum + 2.0d0*fun(x)

        ! Contribution of the remaining points
        else
            sum = sum + 3.0d0*fun(x)
        endif

    enddo

    ! Integral estimation
    integral = sum*(3.0d0*h/8.0d0)

    return

end subroutine simpson_3_8