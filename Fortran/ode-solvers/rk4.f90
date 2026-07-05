!=================================================================
! RUNGE-KUTTA 4 ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Implement the 4th order Runge-Kutta method (RK4) 
!             to solve differential equations
!
! ERROR: O(dt^4) locally, O(dt^3) globally
!
! REQUIREMENTS : External Fortran Subroutine containing the system of
!                diferential equations to be solved is required
!
! INPUTS:
!         ·t    : current time
!         ·dt   : time step
!         ·yyin : input state vector (dimension nequs)
!         ·nequs: number of equations (dimension of the system)
!
! OUTPUTS:
!         ·y_out: output state vector after time dt
!=================================================================

SUBROUTINE MIRUNGEKUTTA4(t, dt, yyin, nequs, y_out)

    IMPLICIT NONE

    ! Inputs
    integer, intent(in) :: nequs
    double precision, intent(in) :: t, dt
    double precision, intent(in) :: yyin(nequs)

    ! Outputs
    double precision, intent(out) :: y_out(nequs)

    ! Internal Variables
    double precision :: yaux(nequs)
    double precision :: K(nequs,4)
    integer :: i

    ! Calculate the K_i of the RK4 algorithm using an auxiliary vector
    ! Requires an external subroutine containing the system of differential equations

    ! K1
    call derivad(t, yyin, K(:,1), nequs)

    ! K2
    do i=1,nequs
        yaux(i) = yyin(i) + 0.5d0*dt*K(i,1)
    end do
    call derivad(t+0.5d0*dt, yaux, K(:,2), nequs)

    ! K3
    do i=1,nequs
        yaux(i) = yyin(i) + 0.5d0*dt*K(i,2)
    end do
    call derivad(t+0.5d0*dt, yaux, K(:,3), nequs)

    ! K4
    do i=1,nequs
        yaux(i) = yyin(i) + dt*K(i,3)
    end do
    call derivad(t+dt, yaux, K(:,4), nequs)

    ! Calculate y_out using the RK4 algorithm
    do i=1,nequs
        y_out(i) = yyin(i) + (dt*(K(i,1)+2*K(i,2)+2*K(i,3)+K(i,4)))/6.0d0
    end do

    RETURN

END SUBROUTINE MIRUNGEKUTTA4