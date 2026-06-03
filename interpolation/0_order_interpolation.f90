!=================================================================
! 0 ORDER INTERPOLATION
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/tuusuario/tu-repo
!=================================================================
! OBJECTIVES: Interpolating f(x) at order 0 from a set of points 
!             {(xi,fi)}
!                       
! INPUTS:
!         ·X,Fi : two vectors containing the points in a matching increasing order
!         ·x    : the point where f(x) is to be interpolated
!
! OUTPUTS: 
!         ·fx   : 0 order interpolation of f at point x
!=================================================================

SUBROUTINE INTERPOLATION_0(x, fx, Xi, Fi)

    implicit none

    ! Inputs
    double precision, intent(in) :: x
    double precision, intent(in) :: Xi(:), Fi(:)

    ! Outputs
    double precision, intent(out) :: fx

    ! Internal variables
    integer :: i, npts

    ! Check both vectors have the same dimension
    if (size(Xi) .eq. size(Fi)) then

        ! Iterate through Xi
        npts = size(xi)
        do i = 1, npts-1

            ! Find within which points x lies
            if (x .ge. Xi(i) .and. x .le. Xi(i+1)) then

                ! Interpolation at order 0
                fx = Fi(i)

                exit

            endif

        enddo

    ! Error when input vectors dimension does not match
    else

        write(*,*) "ERROR: The dimension of the input arrays does not match"

    endif

    RETURN

END SUBROUTINE INTERPOLATION_0
