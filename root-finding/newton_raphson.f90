!=================================================================
! NEWTON-RAPHSON METHOD
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Finding a root of f(x) = 0 using the Newton-Raphson
!             algorithm
!
! REQUIREMENTS: f'(x0) not being too close to zero
!
! CONVERGENCE ORDER: Quadratic convergence near the root
!                       
! INPUTS:
!         ·x0, the point where the iteration begins
!         ·f(x) defined as an external function [fun(x,f(x),f'(x)]
!         ·Maximum number of iterations permited
!         ·Error bound
!
! OUTPUTS: 
!         ·Estimated root within the error range
!         ·Number of iterations needed
!=================================================================

SUBROUTINE NEWTONRAP(x0,er_bound,fun,niter,root,maxiter)

    implicit none

    ! Inputs
    double precision, intent(in) :: x0, er_bound
    integer, intent(in) :: maxiter
    external :: fun 

    ! Output
    double precision, intent(out) :: root
    integer, intent(out) :: niter

    ! Internal variables
    double precision :: f, df, x, dif
    integer :: i

    ! Inicialize the variables
    niter = 0
    x = x0
    dif = 1.0d6 ! This variable controls the error convergence, so it is inicialized in a big number to allow the algorithm to begin 

    ! The method can diverge, so a maximum iteration number will be defined in order to stop it
    do i=1,maxiter

        ! Calculate f(x) and f'(x)
        call fun(x,f,df)

        ! NEWTON-RAPHSON ALGORITHM
        root = x - (f/df)
        dif = abs(root - x)

        ! Stop the program once within the wished error bound
        if (dif < er_bound) then

            niter = i
            x = root
            EXIT

        ! If f'(x) is too small the alorithm will diverge and we must stop the iterations
        elseif (abs(df) < 1.0d-13) then

            write(*,*) "DIVERGENCE ERROR: f'(x) too small"
            EXIT

        else

            ! Prepare the next iteration
            x = root
            niter = i

        endif 

    enddo

    ! Error message when the maximum number of iterations is passed without convergence
    if (i > maxiter .and. dif > er_bound) then
        write(*,*) "ERROR: Maximum number of iterations reached"
    endif

    RETURN

END SUBROUTINE NEWTONRAP
