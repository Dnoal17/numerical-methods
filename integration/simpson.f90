!=================================================================
! SIMPON ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Integrating f(x) in the interval [a,b] through the 
!             Simpson algorithm
!
! ERROR: O(h^4), with h the length of the subintervals 
!                     
! INPUTS:
!         ·a,b the points that define the integration interval
!         ·k an integer that will define the number of subintervals used (N=2^k)
!         ·f(x) defined as an external function [fun(x)]
!
! OUTPUTS: 
!         ·Integral aproximation
!=================================================================

SUBROUTINE SIMPSON(a,b,k,fun,integral)

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

    ! Inicialize the variables
    sum = 0.0d0
    N = 2**k
    h = (b-a)/dble(N)

    ! Iterate through the number of subintervals
    do i=0,N

        ! Define the points
        x = a + h*dble(i)

        ! SIMSON ALGORITHM

        !Contribution of the endpoints
        if (i == 0 .or. i == N) then
            sum = sum + fun(x) 

        !Contribution of even-labled points
        else if (mod(i,2) == 0) then
            sum = sum + 2.0d0*fun(x)

        !Contribution of odd-labled points
        else
            sum = sum + 4.0d0*fun(x)
        endif

    enddo

    ! Integral estimation
    integral = sum*(h/3.0d0)

    RETURN

END SUBROUTINE SIMPSON
