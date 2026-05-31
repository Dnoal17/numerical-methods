!--------------
!NEWTON RAPHSON
!--------------

!Aquesta subrutina implementa el mètode de Newton Raphson per trobar arrels de funcions
!Aquest mètode és molt més eficient que la bisecció, però s'ha de tenir cura d'implementar-lo en punts amb derivades suficientment diferents de 0
!Si la recta tangent és molt horitzontal (df(x0) ~ 0) el mètode no convergeix
SUBROUTINE NEWTONRAP(x0,eps,fun,niter,xarell)

    IMPLICIT NONE

    double precision :: x0, eps, xarell, f, df, x, dif 
    integer :: niter, maxiter = 1000, i
    external :: fun !Funció a trobar les arrels, definida com a Subroutine

    !Incialitzem el comptador d'iteracions i el punt on començarem el mètode
    niter = 0
    x = x0

    !Inicialitzem en un la distància entre l'aproximació de l'arrel i el punt on calculàvem la darrera iteració
    !Ho fem en un valor absurdament gran perque la primera iteració no pugui donar-se per bona mai
    dif = 1.0e6 

    !Com el mètode pot no convergir no utilitzarem un while, si no un do amb un nombre màxim d'iteracions raonable
    do i=1,maxiter

        !Calculem les imatges i derivades de la funció en el punt iterat
        call fun(x,f,df)

        !Si la derivada és molt petita (menos que la preccisió de fortran) sortim del mètode i avisem per pantalla
        if (abs(df) < 1.0d-13) then
            write(*,*) "ERROR: El punt inicial té una derivada molt petita"
            
            EXIT

        !En assolir la precissió dessitjada sortim del bucle amb el darrer valor de Xarell
        elseif (dif < eps) then

            EXIT

        else

            !Algoritme de Newton Raphson
            xarell = x - (f/df)
            dif = abs(x-xarell)

            !Preparem la següent iteració
            x = xarell
            niter = niter + 1

        endif 

    enddo

    !Si el mètode no convergeix en el nombre màxim d'iteracions notifiquem per pantalla
    if (i == maxiter .and. dif > eps) then
        write(*,*) "S'ha arribat al nombre màxim d'iteracions permeses, no s'ha assolit a la preccisió dessitjada"
    endif

    RETURN

END SUBROUTINE NEWTONRAP
