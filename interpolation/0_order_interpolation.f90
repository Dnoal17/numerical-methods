
!----------------------
!INTERPOLACIÓ D'ORDRE 0
!----------------------

!Aquesta subrutina retorna el valor interpolat d'ordre 0 d'un punt x(t) donada un conjunt de fi(xi)
SUBROUTINE XINTERPO0(t, x, npts)

    IMPLICIT NONE

    integer :: i, npts
    real :: t, x
    real :: TI(npts), XI(npts)

    !S'espera que els vectors TI i XI tinguin el conjunt de temps i posicions de la funció a interpolar (podrian posar-se com arguments)
    COMMON /DADES/ TI, XI

    !La interpolació és d'ordre 0, és a dir, es dona al punt interpolat el valor del punt inmediatament anterior
    !Calculem entre quins punt està i hi assignem el valor corresponent

    do i = 1, npts-1
        if (t .ge. TI(i) .and. t .le. TI(i+1)) then

            !Interpolació
            x = XI(i)

        endif
    enddo

    RETURN

END SUBROUTINE XINTERPO0
