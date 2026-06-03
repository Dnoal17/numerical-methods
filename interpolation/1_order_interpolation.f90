!=================================================================
! LINEAR INTERPOLATION
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Interpolating f(x) at first order from a set of points 
!             {(xi,fi)}
!                       
! INPUTS:
!         ·X,Fi : two vectors containing the points in a matching increasing order
!         ·x    : the point where f(x) is to be interpolated
!
! OUTPUTS: 
!         ·fx   : 1st order interpolation of f at point x
!=================================================================

SUBROUTINE LINEAR_INTERPOLATION(x, fx, Xi, Fi)

    implicit none

    ! Inputs
    double precision, intent(in) :: x
    double precision, intent(in) :: Xi(:), Fi(:)

    ! Outputs
    double precision, intent(out) :: fx

    ! Internal variables
    integer :: i, npts
    double precision :: m, n

    ! Check both vectors have the same dimension
    if (size(Xi) .eq. size(Fi)) then

        ! Iterate through Xi
        npts = size(Xi)
        do i = 1, npts-1

            ! Find within which points x lies
            if (x .ge. Xi(i) .and. x .le. Xi(i+1)) then

                ! Linear interpolation (slope and intercept)
                m = (Fi(i+1) - Fi(i)) / (Xi(i+1) - Xi(i))
                n = Fi(i) - m*Xi(i)

                ! Interpolated value
                fx = m*x + n

                exit

            endif

        enddo

    else

        write(*,*) "ERROR: The dimension of the input arrays does not match"

    endif

    RETURN

END SUBROUTINE LINEAR_INTERPOLATION