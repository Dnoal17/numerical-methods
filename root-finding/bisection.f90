!---------------------
!MÈTODE DE LA BISECCIÓ
!---------------------

!Aquesta subrutina implementa el mètode de la bisecció per trobar arrels de funcions
!Aquest, és un mètode comptacte: Troba una arrels sempre, encara que no sigui el més eficient
SUBROUTINE BISECTION(A,B,eps,fun,niter,xarell)

    IMPLICIT NONE

    double precision :: A, B, C, eps, xarell
    double precision :: fA, dfA, fB, dfB, fC, dfC
    integer :: niter, i
    external :: fun !Funció a trobar les arrels, definida com a Subroutine

    !Incialitzem els valors de fA,fB per iniciar l'algoritme
    !Les variables dfA, dfB no fan res aquí, pero a la pràctica era convenient definir les funcions amb les derivades per reutilizar-les en altres mètodes
    call fun(A,fA,dfA)
    call fun(B,fB,dfB)

    !El mètode de la bisecció permet saber el nombre de passos necessàris per trobar una arrel amb precissió eps
    !Això és possible perquè coneixem el tamany de l'interval en la n-éssima iteració (es redueix 1/2 amb cadascuna)
    niter = int(dlog((B-A)/eps)/dlog(2.0d0)) + 1

    !Inciem l'algoritme comprobant que hi ha mínim una arrel en l'intèrval seleccionat (Th. Bolzano)
    if (fA*fB<0) then

        do i=1,niter
            
            !Definim C i calculem la seva imatge
            C = (A+B)/2
            call fun(C,fC,dfC)

            !Continuem el mètode comprovant en quina partició de l'intèrval roman l'arrel
            if (fC*fA<0.0) then

                B = C
                fB = fC

            elseif (fC*fB<0.0) then

                A = C
                fA = fC

            endif

        enddo

        !En acabar el mètode, agafem com arrel l'ultim valor de C
        xarell = C

    !Si hem seleccionat malament els valors de fA i fB ho escribim per pantalla per assabentar-nos
    elseif (fA*fB>0) then
        write(*,*) "NO HI HA CANVIA DE SIGNE ENTRE A =", A, "I B =", B

    !Potser encertem l'arrel en donar els extrems de l'intèrval i ja hem acabat (poc probable)
    elseif (fA == 0.0) then
            C = A
            xarell = C
    elseif (fB ==0.0) then
            C = B
            xarell = C
    endif

    RETURN

END SUBROUTINE BISECTION