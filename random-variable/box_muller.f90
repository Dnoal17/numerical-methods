!=================================================================
! BOX-MULLER METHOD
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Generate Gaussian random numbers N(0,1) using the 
!             Box-Müller transformation method
!
! INPUTS:
!         ·ndat  : number of data points to generate
!
! OUTPUTS:
!         ·xgaus : array containing the generated Gaussian numbers N(0,1)
!
! NOTES:
!         To generate Gaussian numbers with other constants sigma, mu:
!         N(mu,sigma) = sigma*N(0,1) + mu
!=================================================================

SUBROUTINE SUBGAUSS(ndat,xgaus)

    implicit none

    ! Inputs
    integer, intent(in) :: ndat

    ! Outputs
    double precision, intent(out) :: xgaus(ndat)

    ! Internal Variables
    double precision :: X1, X2, r, theta
    double precision, parameter :: Pi = dacos(-1.0d0)
    integer :: i

    ! BOx-MÜLLER METHOD
    do i=1,ndat,2

        ! Generate 2 random numbers in [0,1]
        X1 = rand() 
        X2 = rand()

        ! Use the variable change of the method
        r = dsqrt(-2.0d0*log(X1))
        theta = 2.0d0*Pi*X2

        ! Obtain two normal numbers N(0,1)
        xgaus(i) = r*dcos(theta)

        ! The second normal number is only stored if possible (if ndat is odd, the last one will be left over)
        if (i+1 .le. ndat) then 
            
            xgaus(i+1) = r*dsin(theta)

        endif

    enddo

    RETURN

END SUBROUTINE SUBGAUSS