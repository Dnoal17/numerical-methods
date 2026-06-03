!=================================================================
! BISECTION METHOD
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Finding a root of f(x) = 0 in the interval [A,B] 
!             using the bisection method
!
! REQUIREMENTS: f(A)·f(B) < 0
!
! CONVERGENCE ORDER: The error in this method is bounded by the
!                    number of iterations as: error <= (B-A)/2^n
!                       
! INPUTS:
!         ·A,B
!         ·f(x) defined as an external function [fun(x,f(x)]
!         ·Error bound
!
! OUTPUTS: 
!         ·Estimated root within the error range
!         ·Number of iterations needed
!=================================================================

SUBROUTINE BISECTION(A,B,er_bound,fun,niter,root)

    implicit none

    ! Inputs
    double precision, intent(in) :: A, B, er_bound
    external :: fun 

    ! Output
    double precision, intent(out) :: root
    integer, intent(out) :: niter

    ! Internal variables
    double precision :: C
    double precision :: fA, fB, fC
    integer :: i

    ! Inicialize f(B) and f(A)
    call fun(A,fA)
    call fun(B,fB)

    ! Number of iterations needed to achieve the error bound
    niter = int(dlog((B-A)/er_bound)/dlog(2.0d0)) + 1

    ! Check the sign requirement is fulfilled
    if (fA*fB<0) then

        ! BISECTION ALORITHM

        ! Iterate for the previously calculated amount
        do i=1,niter
            
            ! Define the midpoint of the interval and calculate its image
            C = (A+B)/2
            call fun(C,fC)

            ! Find which new interval, [A,C] or [C,B], is the root in by checking where there is a change of sign

            if (fC*fA<0.0) then

                ! Prepare the next iterations
                B = C
                fB = fC

            elseif (fC*fB<0.0) then

                ! Prepare the next iterations
                A = C
                fA = fC

            endif

        enddo

        ! Take the last midpoint as the root estimation
        root = C

    ! If the input does not fulfill the requirement notify the reasson why 
    elseif (fA*fB>0) then
        write(*,*) "ERROR: Both A =", A, " and B =", B "have the same sign"

    ! It may be the case where A or B are the root themselves
    elseif (fA == 0.0) then
            C = A
            root = C
    elseif (fB ==0.0) then
            C = B
            root = C
    endif

    RETURN

END SUBROUTINE BISECTION