!----------------------------
!MÈTODE D'ACCEPTACIÓ I REBUIG
!----------------------------

!Aquesta subrutina implementa el mètode d'acceptació i rebuig per aconseguir números aleatoris distribuits segons una densitat de probabilitat g(x)
!És important que la densitat de probabilitat estigui definida (i normalitzada) en els mateixos a,b en els que generarem els números
SUBROUTINE SUBAIR(ndat, xnums, fun, a, b, M)

    IMPLICIT NONE

    integer :: ndat, i, iseed
    double precision :: M, a, b, x, p, xnums(ndat)
    double precision, external :: fun !Densitat de probabilitat

    !Calculem la seed amb el meu NIUB
    iseed = 20795946
    call srand(iseed)

    !Inicialitzem l'output
    xnums = 0.0d0

    !Generem els nombres aleatoris donada una densitat de probabilitat f(x)
    i = 1

    !Utilitzem un bucle que no s'acaba fins tenir ndat números
    do while (i .le. ndat)

    !Es generen 2 números aleatoris entre [0,1]
        x = rand() 
        p = rand()

        !Convertim aquests números en uniformes entre [a,b] i [0,M]
        x = a + x * (b - a)
        p = p * M

        !Condició d'acceptació del mètode
        if (fun(x) .ge. p) then

            xnums(i) = x

            !Si acceptem un número contem +1 al bucle
            i = i + 1

        endif

    enddo

    RETURN

END SUBROUTINE SUBAIR