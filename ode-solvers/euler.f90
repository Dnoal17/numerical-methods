
!--------------------
!MÈTODE D'EULER BÀSIC
!--------------------

!Subrutina que aplica el mètode d'Euler per resoldre equacions diferencials plantejant-les com un sistema d'equacions
!Requeries funcions externes que continguin la informació del sistema d'equacions
SUBROUTINE METODE_EULER(N, y, fun, a, b, nvar)

    IMPLICIT NONE

    integer :: i, j, N, nvar
    double precision :: a, b, h, t
    double precision :: y(N, nvar)
    double precision, external :: fun !Fun conté un sistema d'equacions amb les derivades de la EDO

    !Determinem el pas h
    h = (b - a)/dble(N - 1)

    !Utilitzem l'algoritme del mètode de Euler. Comencem el bucle en i=2 perquè partim de les condicions inicials
    do i = 2, N

        !Calculem el valor de la variable independent per cada iteració
        t = a + h*dble(i-2)

        !Actualitzem cada variable del sistema d'equacions diferencials. 
        !La variable j té la funció de canviar d'equació dins el sistema d'equacions, permetent resoldre equacions de segon ordre (que són les que ens interesen)
        do j = 1, nvar

            y(i,j) = y(i-1,j) + h * fun(t, y(i-1,:), j)

        end do

    end do

    RETURN

END SUBROUTINE METODE_EULER