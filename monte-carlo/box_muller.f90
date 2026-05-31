!-----------------
!MÈTODE BOX-MÜLLER
!-----------------

!Aquesta subrutina implementa el mètode Box-Müller per generar números aleatoris gaussians N(0,1)
!Per generar números gaussians amb altres constants sigma,mu utilitzem que N(mu,sigma) = sigma*N(0,1) + mu
SUBROUTINE SUBGAUSS(ndat,xgaus)

    implicit none
    double precision :: xgaus(ndat), X1, X2, r, theta 
    double precision, parameter :: Pi = dacos(-1.0d0)
    integer :: ndat, i, iseed = 20795946

    !Calculem la seed amb el meu NIUB
    call srand(iseed)

    !Apliquem el mètode Box-Müller
    do i=1,ndat,2

        !Generem 2 números aleatoris en [0,1]
        X1 = rand() 
        X2 = rand()

        !Utilitzem el canvi de variable del mètode
        r = dsqrt(-2.0d0*log(X1))
        theta = 2.0d0*Pi*X2

        !Trobem dos números normals N(0,1)
        xgaus(i) = r*dcos(theta)

        if (i+1 .le. ndat) then 
            
            !El segon número normal el posem només si podem, ja que si ndat és senar sobrarà l'ultim
            xgaus(i+1) = r*dsin(theta)

        endif

	enddo

    RETURN

END SUBROUTINE SUBGAUSS
