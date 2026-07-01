!=================================================================
! SECANT METHOD
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Finding a root of f(x) = 0 using the secant method.
!
! REQUIREMENTS: abs(f(x1)-f(x0)) being not to close to zero
!
! CONVERGENCE ORDER: Superlineal convergence near the root
!                       
! INPUTS:
!         ·x0,x1    : points where the iteration begins
!         ·fun      : defined as an external function [fun(x,f(x))]
!         ·mwxiter  : maximum number of iterations permited
!         ·er_bound : desired error bound
!
! OUTPUTS: 
!         ·root     : estimated root within the error range
!         ·niter    : number of iterations needed
!=================================================================


SUBROUTINE SECANT(x0,x1,er_bound,fun,niter,root,maxiter)

    implicit none
    
    ! Inputs
    double precision, intent(in) :: x0, x1, er_bound
    integer, intent(in) :: maxiter
    external :: fun 

    ! Output
    double precision, intent(out) :: root
    integer, intent(out) :: niter

    ! Internal variables
    double precision :: f0, f1, x, xprev, dif
    integer :: i

    ! Inicialize the variables
    niter = 0
    xprev = x0
    x = x1
    root = x1
    dif = 1.0d6 ! This variable controls the error convergence, so it is inicialized in a big number to allow the algorithm to begin 

    ! The method can diverge, so a maximum iteration number will be defined in order to stop it
    do i = 1, maxiter

        ! Calculate f(x1) and f(x0) each iteration
        call fun(xprev,f0)
        call fun(x,f1)

        ! Check the method will not diverge
        if (abs(f1 - f0) < 1.0d-13) then

            write(*,*) "ERROR: Divergence error "

            EXIT

        endif

        ! SECANT ALGORITHM
        root = x - f1 * (x - xprev) / (f1 - f0)
        dif = abs(x - root)

        ! Prepare the next iteration
        xprev = x
        x = root
        niter = niter + 1

        ! Check if the error bound has been achieved
        if (dif < er_bound) then

            EXIT

        endif

    enddo

    ! If the method could not converge within the maximum number of iterations send an error message
    if (i == maxiter .and. dif > er_bound) then
        write(*,*) "ERROR: Maximum number of iterations was reached without achieving the error bound"
    endif

    RETURN

END SUBROUTINE SECANT
