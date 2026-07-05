!=================================================================
! MOMENT CALCULATOR
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Calculate the m-th order central moment, expected value 
!             and standard deviation of a set of random numbers 
!
! INPUTS:
!         ·ndat           : number of data points
!         ·xdata          : array containing the data
!         ·m              : order of the central moment to calculate
!
! OUTPUTS:
!         ·expected_x     : expected value
!         ·moment_m       : m-th order central moment
!         ·stnd_deviation : standard deviation (only if m=2, otherwise 0)
!=================================================================

subroutine moment_order_m(ndat, xdata, m, expected_x, moment_m, stnd_deviation)
    
    implicit none

    ! Inputs
    integer, intent(in) :: ndat, m
    double precision, intent(in) :: xdata(ndat)

    ! Outputs
    double precision, intent(out) :: expected_x, moment_m, stnd_deviation

    ! Internal Variables
    integer :: i

    ! Initialize the variables
    expected_x = 0.0d0
    moment_m = 0.0d0

    ! Calculate the expected value
    do i=1,ndat

        expected_x = xdata(i) + expected_x

    enddo

    expected_x = expected_x/ndat

    ! Calculate the m-th order moment
    do i=1,ndat

        moment_m = (xdata(i) - expected_x)**m + moment_m

    enddo

    moment_m = moment_m/ndat

    ! If the moment is second order, also calculate the standard deviation
    if (m .eq. 2) then
        stnd_deviation = sqrt(moment_m)
    else
        stnd_deviation = 0.0d0
    endif

    return

end subroutine MOMENT_ORDER_M