!=================================================================
! MONTECARLO ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Integrating f(x) in N-dimensional interval [a,b] 
!             through the Monte Carlo algorithm
!
! ERROR: O(1/√N), with N the number of random points
!                     
! INPUTS:
!         ·a,b vectors that define the integration interval in N dimensions
!         ·dim integer that defines the number of dimensions
!         ·N an integer that defines the number of random points used
!         ·f(x) defined as an external function [fun(x)]
!
! OUTPUTS: 
!         ·Integral aproximation
!         ·Error estimation
!=================================================================

SUBROUTINE CRUDE_MONTECARLO(dim,N,fun,a,b,integral,error)

    implicit none

    ! Inputs
    integer, intent(in) :: dim, N
    double precision, intent(in) :: a(dim), b(dim)
    external :: fun

    ! Outputs
    double precision, intent(out) :: integral, error

    ! Internal Variables
    double precision :: lambda, expected_x, moment_2, sigma
    double precision :: x(N,dim), f(N)
    integer :: k, j

    ! Initialize the variables
    lambda = 1.0d0
    integral = 0.0d0

    ! Generate N*dim random numbers
    do k=1,N

        ! Fill each row with dim random numbers uniformly distributed
        do j=1,dim
            x(k,j) = rand()
            x(k,j) = a(j) + (b(j) - a(j))*x(k,j)
        enddo

        ! Save the images of each point
        f(k) = fun(x(k,:))

        ! Sum all contributions
        integral = integral + f(k)

    enddo

    ! Determine the (b-a) factor from the method theory
    do j=1,dim
        lambda = lambda*(b(j)-a(j))
    enddo

    ! Monte Carlo approximation of the integral
    integral = (lambda/N)*integral

    ! Calculate expectation and second moment
    expected_x = 0.0d0
    moment_2  = 0.0d0

    do k=1,N
        expected_x = expected_x + f(k)
        moment_2  = moment_2  + f(k)**2
    enddo

    expected_x = expected_x / N
    moment_2 = moment_2 / N

    ! Standard deviation of the random variable f(x)
    sigma = dsqrt(moment_2 - expected_x**2)

    ! Error of the Monte Carlo method: σ/√N multiplied by the domain volume
    error = (lambda * sigma) / dsqrt(dble(N))

END SUBROUTINE CRUDE_MONTECARLO