!-------------------
!MÈTODE DE LA SECANT
!-------------------

!Aquesta subrutina implementa el mètode de la secant per trobar arrels de funcions
!Aquest mètode és més eficient que la bisecció i menys que Newton Raphson, però s'ha de tenir cura d'implementar-lo en punts amb derivades suficientment diferents de 0
!Si la recta tangent és molt horitzontal (df(x0) ~ 0) el mètode no convergeix
SUBROUTINE SECANT(x0,x1,eps,fun,niter,xarell)

    IMPLICIT NONE

    double precision :: x0, x1, eps, xarell
    double precision :: f0, f1, df0, df1
    double precision :: x, xprev, dif
    integer :: niter, maxiter = 1000, i
    external :: fun   !Funció a trobar les arrels, definida com a Subroutine

    !Inicialitzem el comptador d'iteracions
    niter = 0

    !Inicialitzem els dos primers punts del mètode
    xprev = x0
    x = x1

    !Inicialitzem la diferència entre iteracions amb un valor gran
    dif = 1.0d6

    !Bucle principal amb nombre màxim d'iteracions
    do i = 1, maxiter

        !Calculem les imatges de la funció als dos punts
        call fun(xprev, f0, df0)
        call fun(x, f1, df1)

        !Comprovem que el denominador no sigui massa petit
        if (abs(f1 - f0) < 1.0d-13) then
            write(*,*) "ERROR: Denominador molt petit en el mètode de la secant"
            EXIT

        !Comprovem si ja hem assolit la precisió desitjada
        elseif (dif < eps) then

            EXIT

        else

            !Algoritme del mètode de la secant
            xarell = x - f1 * (x - xprev) / (f1 - f0)

            !Calculem la diferència entre iteracions consecutives
            dif = abs(x - xarell)

            !Preparem la següent iteració
            xprev = x
            x = xarell
            niter = niter + 1

        endif

    enddo

    !Si el mètode no convergeix en el nombre màxim d'iteracions
    if (i == maxiter .and. dif > eps) then
        write(*,*) "S'ha arribat al nombre màxim d'iteracions permeses, no s'ha assolit la precisió desitjada"
    endif

    RETURN

END SUBROUTINE SECANT
