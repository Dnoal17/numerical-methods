!=================================================================
! ACCEPTANCE-REJECTION METHOD
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Generate random numbers following a given probability 
!             density function g(x) using the acceptance-rejection method
!
! INPUTS:
!         ·ndat  : number of data points to generate
!         ·fun   : probability density function (external function)
!         ·a     : lower bound of the interval [a,b]
!         ·b     : upper bound of the interval [a,b]
!         ·M     : upper bound such that fun(x) <= M for all x in [a,b]
!
! OUTPUTS:
!         ·xnums : array containing the generated random numbers
!
! NOTES:
!         The probability density function must be defined (and normalized)
!         in the same interval [a,b] where the numbers will be generated
!=================================================================

SUBROUTINE SUBAIR(ndat, xnums, fun, a, b, M)

    IMPLICIT NONE

    integer :: ndat, i
    double precision :: M, a, b, x, p, xnums(ndat)
    double precision, external :: fun

    ! Initialize output array
    xnums = 0.0d0

    ! Generate random numbers following the probability density function f(x)
    i = 1

    ! Use a loop that doesn't stop until we have ndat numbers
    do while (i .le. ndat)

        ! Generate 2 random numbers in [0,1]
        x = rand() 
        p = rand()

        ! Convert these numbers to uniform distributions in [a,b] and [0,M]
        x = a + x * (b - a)
        p = p * M

        ! Acceptance condition of the method
        if (fun(x) .ge. p) then

            xnums(i) = x

            ! If we accept a number, count +1 in the loop
            i = i + 1

        endif

    enddo

    RETURN

END SUBROUTINE SUBAIR