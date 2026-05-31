! Method: Linear interpolation
! Input: t, time and position arrays
! Output: interpolated value x(t)

!-------------------
!INTERPOLACIÓ LINEAL
!-------------------

!Aquesta subrutina retorna el valor interpolat linealment d'un punt x(t) donada un conjunt de fi(xi)
SUBROUTINE XINTERPO(t, x, npts)

    IMPLICIT NONE

    integer :: i, npts
    real :: t, x
    real :: TI(npts), XI(npts)

    real :: m, n

    !S'espera que els vectors TI i XI tinguin el conjunt de temps i posicions de la funció a interpolar (podrian posar-se com arguments)
    COMMON /DADES/ TI, XI

    !La interpolació és lineal, és a dir, es situa el punt interpolat sobre la recta que uneix els dos punts més propers
    !Calculem quins són aquest punts i fem la recta
    do i = 1, npts-1
        if (t .ge. TI(i) .and. t .le. TI(i+1)) then

            m = (XI(i+1) - XI(i)) / (TI(i+1) - TI(i))
            n = XI(i) - m*TI(i)

            !Interpolació
            x = m*t + n

        endif
    enddo

    RETURN

END SUBROUTINE XINTERPO