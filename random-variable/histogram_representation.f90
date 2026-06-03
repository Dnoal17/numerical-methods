!=================================================================
! HISTOGRAM REPRESENTATION
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Constructing a normalized histogram from a data set
!             and estimating the statistical uncertainty of each bin
!
! INPUTS:
!         ·ndat  : number of data points
!         ·xdata : array containing the data sample
!         ·nbox  : number of histogram bins
!
! OUTPUTS:
!         ·xhisto   : center of each histogram bin
!         ·histo    : normalized histogram (probability density)
!         ·errhisto : statistical uncertainty of each bin
!=================================================================

SUBROUTINE HISTOGRAM(ndat,xdata,nbox,xhisto,histo,errhisto)

    implicit none

    ! Inputs
    integer, intent(in) :: nbox, ndat
    double precision, intent(in) :: xdata(ndat)

    ! Outputs
    double precision, intent(out) :: xhisto(nbox), histo(nbox), errhisto(nbox)

    ! Internal Variables
    integer :: i_caixa, i
    double precision :: x_min, x_max, h, p_i
    
    ! Determine the minimum and maximum values of the data sample
    x_min = minval(xdata)
    x_max = maxval(xdata)

    ! Define the bin width using the selected number of bins
    h = (x_max-x_min)/nbox

    ! Initialize all output arrays
    histo = 0.0d0
    xhisto = 0.0d0
    errhisto = 0.0d0

    ! Identify the bin corresponding to each data point
    do i=1,ndat

        ! The maximum value must be treated separately to avoid
        ! assigning it to a non-existent bin
        if (xdata(i) .eq. x_max) then
            i_caixa = nbox

        ! General case
        else
            i_caixa = int((xdata(i)-x_min)/h) + 1

        endif

        ! Increase the counter of the corresponding bin
        histo(i_caixa) = histo(i_caixa) + 1

    enddo

    do i=1,nbox

        ! Compute the center of each bin
        xhisto(i) = x_min + (i-0.5d0)*h

        ! Normalize the histogram to obtain a probability density
        histo(i) = histo(i)/(ndat*h)

        ! Estimate the statistical uncertainty assuming binomial statistics in each bin
        p_i = histo(i)*h
        errhisto(i) = dsqrt(p_i*(1.0d0-p_i)/ndat)/h

    enddo

    RETURN

END SUBROUTINE HISTOGRAMA