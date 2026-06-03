!=================================================================
! POISSON EQUATION SOLVER ALGORITHM
!=================================================================
! Author: Daniel Noal Pineda
! Email : noaldaniel41@gmail.com
! Date  : 2025
! Repository: https://github.com/Dnoal17/numerical-methods.git
!=================================================================
! OBJECTIVES: Solve the Heat Equation (2D Poisson equation) using Jacobi, 
!             Gauss-Seidel, and Successive Over-Relaxation (SOR) methods
!
! NOTES: The method used is controlled by the variable icontrol
!        For correct operation, must be implemented within a CONTAINS 
!        block (does not work properly at the end of the file)
!
! INPUTS:
!         ·icontrol     : method selector (1=Jacobi, 2=Gauss-Seidel, 3=SOR)
!         ·t_cinitials  : initial temperature matrix
!         ·p            : source term matrix
!         ·Lx           : domain length in x-direction
!         ·Ly           : domain length in y-direction
!         ·h            : grid spacing
!         ·i_point      : i-index for convergence monitoring point
!         ·j_point      : j-index for convergence monitoring point
!         ·omega        : relaxation parameter (for SOR method)
!         ·unit_file    : file unit for writing convergence data
!
! OUTPUTS:
!         ·T_final      : final temperature matrix (solution)
!         ·conv_p       : convergence vector at the monitoring point
!=================================================================

SUBROUTINE SOLVER_POISSON_2D(icontrol,t_cinitials,T_final,p,Lx,Ly,h,i_point,j_point,conv_p,omega,unit_file)

    implicit none

    ! Inputs
    integer, intent(in) :: icontrol, i_point, j_point, unit_file
    double precisionn, intent(in) :: t_cinitials(:,:), p(:,:)
    double precisionn, intent(in) :: Lx, Ly, h, omega

    ! Outputs
    double precisionn, intent(out) :: T_final(:,:), conv_p(:)

    ! Internal Variables
    logical :: convergence
    integer :: i, j, Nx, Ny, iteracio, max_iter
    double precisionn :: Error, precision
    double precisionn, allocatable :: T_k(:,:), T_k_plus_1(:,:) 

    ! Convergence parameters
    precision = 1.0d-4
    max_iter = 100000
    
    ! Matrix dimensions
    Nx = size(t_cinitials, 1)
    Ny = size(t_cinitials, 2)

    ! Memory allocation
    allocate(T_k(Nx,Ny), T_k_plus_1(Nx,Ny))
    
    ! Initialization
    convergence = .false.
    iteracio = 0
    
    ! Copy initial conditions
    T_k = t_cinitials
    T_k_plus_1 = T_k

    ! Iterate until convergence
    do while (.not. convergence .and. iteracio < max_iter)
        
        iteracio = iteracio + 1

        !-------------------------------
        ! JACOBI METHOD (ICONTROL = 1)
        !-------------------------------

        if (icontrol == 1) then
            do i = 2, Nx-1
                do j = 2, Ny-1

                    ! Jacobi method algorithm
                    T_k_plus_1(i,j) = (T_k(i+1,j) + T_k(i-1,j) + &
                                        T_k(i,j+1) + T_k(i,j-1) + &
                                        h**2 * p(i,j)) / 4.0d0

                enddo
            enddo

        !-------------------------------------
        ! GAUSS-SEIDEL METHOD (ICONTROL = 2)
        !-------------------------------------

        elseif (icontrol == 2) then
            do i = 2, Nx-1
                do j = 2, Ny-1

                    ! Gauss-Seidel algorithm
                    T_k_plus_1(i,j) = (T_k_plus_1(i-1,j) + T_k(i+1,j) + &
                                      T_k_plus_1(i,j-1) + T_k(i,j+1) + &
                                      h**2 * p(i,j)) / 4.0d0

                enddo
            enddo

        !---------------------------------------
        ! SOR METHOD (ICONTROL = 3)
        !---------------------------------------

        elseif (icontrol == 3) then
            do i = 2, Nx-1
                do j = 2, Ny-1

                    ! Successive Over-Relaxation (SOR) algorithm
                    T_k_plus_1(i,j) = T_k(i,j) + omega * (T_k_plus_1(i-1,j) + T_k(i+1,j) + &
                                      T_k_plus_1(i,j-1) + T_k(i,j+1) - 4*T_k(i,j) + &
                                      h**2 * p(i,j)) / 4.0d0

                enddo 
            enddo
        endif

        !------------------------------------------------
        ! CONVERGENCE - ERROR CALCULATION - NEXT ITERATION
        !------------------------------------------------
        
        ! Study convergence at a point and save data to file
        write(unit_file, '(I20, F20.6)') iteracio, T_k_plus_1(i_point,j_point)

        ! Save to convergence vector
        conv_p(iteracio) = T_k_plus_1(i_point,j_point)

        ! Calculate the error
        Error = 0.0d0

        do i = 2, Nx-1
            do j = 2, Ny-1
                Error = max(Error, abs(T_k_plus_1(i,j) - T_k(i,j)))
            enddo
        enddo

        ! Check convergence
        if (Error < precision) then
            convergence = .TRUE.
        endif

        ! Prepare for the next iteration
        T_k = T_k_plus_1
    
    enddo

    ! Save the final result
    T_final = T_k_plus_1
    
    deallocate(T_k, T_k_plus_1)

    RETURN

END SUBROUTINE SOLVER_POISSON_2D