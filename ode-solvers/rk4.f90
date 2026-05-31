
!-------------
!RUNGE–KUTTA 4
!-------------

!Aquesta subrutina implementa el mètode de RK4 per resoldre equacions diferencials.
!És un mètode molt eficient perque té un error local i, doncs, l'error global no augmenta en augmentar l'interval d'integració
SUBROUTINE MIRUNGEKUTTA4(t, dt, yyin, nequs, yyout)

    IMPLICIT NONE

    integer :: nequs, i
    double precision :: t, dt
    double precision :: yyin(nequs), yyout(nequs), yaux(nequs)
    double precision :: K(nequs,4)

    !Calculem les K_i de l'algoritme RK4 ajudant-nos d'un vector auxiliar
    !Necesitem una subrutina externa que contingui el sitema d'equacions diferencials que es vol resoldre

    !K1
    call derivad(t, yyin, K(:,1), nequs)

    !K2
    do i=1,nequs
        yaux(i) = yyin(i) + 0.5d0*dt*K(i,1)
    end do
    call derivad(t+0.5d0*dt, yaux, K(:,2), nequs)

    !K3
    do i=1,nequs
        yaux(i) = yyin(i) + 0.5d0*dt*K(i,2)
    end do
    call derivad(t+0.5d0*dt, yaux, K(:,3), nequs)

    !K4
    do i=1,nequs
        yaux(i) = yyin(i) + dt*K(i,3)
    end do
    call derivad(t+dt, yaux, K(:,4), nequs)

    !Calculem yyout utitlizant l'algoritme RK4
    do i=1,nequs
        yyout(i) = yyin(i) + (dt*(K(i,1)+2*K(i,2)+2*K(i,3)+K(i,4)))/6.0d0
    end do

    RETURN

END SUBROUTINE MIRUNGEKUTTA4