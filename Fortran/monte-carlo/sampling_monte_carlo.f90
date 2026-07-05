!=================================================================
! MONTECARLO SAMPLING ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Integrating f(x) in N-dimensional interval through 
!             the Monte Carlo algorithm with importance sampling
!             using a probability density function
!
! REQUIREMENT: Density fucntion must be normalized
!
! ERROR: The more similar the function and the probability density 
!        are, the more efficient the method will be.
!
! INPUTS:
!         ·dim     : integral dimension (integeer)
!         ·N       : number of random points used each integration (integeer)
!         ·fun     : function to be integrated defined as an external Fortran function [f(x)]
!         ·density : probability density for sampling, defined as Fortran function [density(x)]
!         ·nums    : matrix containing the random numbers (N x dim)
!
! OUTPUTS: 
!         ·integral : integral approximation
!         ·error    : error estimation
!=================================================================

subroutine sampling_montecarlo(dim, N, fun, density, nums, integral, error)

    implicit none

    ! Inputs
    integer, intent(in) :: dim, N
    double precision, intent(in) :: nums(N,dim)
    external :: fun, density

    ! Outputs
    double precision, intent(out) :: integral, error

    ! Internal Variables
    double precision :: expected_x, moment_2, sigma
    double precision :: f(N), g(N)
    integer :: i

    ! Initialize the variables
    integral = 0.0d0

    ! Save the images of the function and the density
    do i=1,N

        f(i) = fun(nums(i,:))
        g(i) = density(nums(i,:))

        ! Sum contributions according to Monte Carlo with sampling method
        integral = integral + f(i)/g(i)

    enddo

    ! Final determination of the integral
    integral = integral/N

    ! Calculate expectation and second moment
    expected_x = 0.0d0
    moment_2  = 0.0d0

    do i=1,N
        expected_x = expected_x + (f(i)/g(i))
        moment_2  = moment_2  + (f(i)/g(i))**2
    enddo

    expected_x = expected_x / N
    moment_2 = moment_2 / N

    ! Standard deviation of the random variable f(x)
    sigma = dsqrt(moment_2 - expected_x**2)

    ! Error of the Monte Carlo method: σ/√N
    error = sigma / dsqrt(dble(N))

    return

end subroutine sampling_montecarlo