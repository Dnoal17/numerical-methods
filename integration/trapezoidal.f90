!--------------------
!MÈTODE DELS TRAPEZIS
!--------------------

!Aquest és un mètode d'integració numèrica aproximant l'àrea sota la corba d'una funció amb trapezis i calculant les seves àreas
SUBROUTINE METODETRAPEZIS(x1,x2,k,funcio,resul)

    IMPLICIT NONE

    double precision :: x1, x2,resul, h, x, suma
    integer :: k, i, N
    double precision, external :: funcio !Funció a integrar

    !Tamany dels intervals i inicialització de la suma
    N = 2**k !Fem particions de 2**k intervals perque el mètode convergeixi
    h = (x2-x1)/dble(N)
    suma = 0.0d0

    do i=0,N

        x = x1 + h*i

        !El pes dels extrems és la meitat que la resta de punts
        if (i == 0 .or. i == N) then

            suma = suma + funcio(x)/2

        else

            suma = suma + funcio(x)

        endif

    enddo

    resul = suma*h !Obtenim el valor de la integral (Error O(h**2))

    RETURN

END SUBROUTINE METODETRAPEZIS