
    !-----------------------------------------
    !MÈTODES NUMÈRICS PER L'EQUACIÓ DE POISSON
    !-----------------------------------------

    !Aquesta subrutina resol l'equació de poisson en 2d amb els mètodes de Jacobi, Gauss-Seidel i Sobrerelaxació
    !El mètode utitlitzat està controlat per la variable icontrol
    !Per que funcioni correctament s'ha d'implementar dins un CONTAINS, no funciona bé al final del document (no sé per què)

    SUBROUTINE SOLVER_POISSON_2D(icontrol,T_cinicials,T_final,p,punt_i,punt_j,conv_p,omega,unit_file)

        IMPLICIT NONE

        logical :: convergencia
        integer :: icontrol, i, j, Nx, Ny, iteracio, max_iter
        double precision :: Error, precisio
        double precision, intent(in) :: T_cinicials(:,:), p(:,:)
        double precision, intent(in) :: omega
        integer, intent(in) :: punt_i, punt_j, unit_file
        double precision, intent(out) :: T_final(:,:), conv_p(:)

        !Aquestes matrius guarden la iteració k i la k+1 als 3 algoritmes.
        double precision, allocatable :: T_k(:,:), T_k_plus_1(:,:) 

        !-------------------------------------------------------
        !-------------------------------------------------------
        !Paràmetres del problema
        double precision :: Lx, Ly, h
        COMMON/PARAMETRES/ Lx, Ly, h
        !-------------------------------------------------------
        !-------------------------------------------------------

        !Paràmetres de convergència
        precisio = 1.0d-4
        max_iter = 100000
        
        !Dimensions de les matrius
        Nx = size(T_cinicials, 1)
        Ny = size(T_cinicials, 2)

        !Allocació de memòria
        allocate(T_k(Nx,Ny), T_k_plus_1(Nx,Ny))
        
        !Inicialització
        convergencia = .false.
        iteracio = 0
        
        !Copia de les condicions inicials
        T_k = T_cinicials
        T_k_plus_1 = T_k

        !Iterarem fins que el mètode convergeixi
        do while (.not. convergencia .and. iteracio < max_iter)
            
            iteracio = iteracio + 1

            !-------------------------------
            !MÈTODE DE JACOBI (ICONTROL = 1)
            !-------------------------------

            if (icontrol == 1) then
                do i = 2, Nx-1
                    do j = 2, Ny-1

                        !Algoritme del mètode de Jacobi
                        T_k_plus_1(i,j) = (T_k(i+1,j) + T_k(i-1,j) + &
                                            T_k(i,j+1) + T_k(i,j-1) + &
                                            h**2 * p(i,j)) / 4.0d0

                    enddo
                enddo

            !-------------------------------------
            !MÈTODE DE GAUSS-SEIDEL (ICONTROL = 2)
            !-------------------------------------

            elseif (icontrol == 2) then
                do i = 2, Nx-1
                    do j = 2, Ny-1

                        !Algoritme de Gauss-Seidel
                        T_k_plus_1(i,j) = (T_k_plus_1(i-1,j) + T_k(i+1,j) + &
                                          T_k_plus_1(i,j-1) + T_k(i,j+1) + &
                                          h**2 * p(i,j)) / 4.0d0

                    enddo
                enddo

            !---------------------------------------
            !MÈTODE DE SOBRERELAXACIÓ (ICONTROL = 3)
            !---------------------------------------

            elseif (icontrol == 3) then
                do i = 2, Nx-1
                    do j = 2, Ny-1

                        !Algoritme de sobrerelaxació (SOR)
                        T_k_plus_1(i,j) = T_k(i,j) + omega * (T_k_plus_1(i-1,j) + T_k(i+1,j) + &
                                          T_k_plus_1(i,j-1) + T_k(i,j+1) - 4*T_k(i,j) + &
                                          h**2 * p(i,j)) / 4.0d0


                    enddo 
                enddo
            endif

            !------------------------------------------------
            !CONVERGÈNCIA - CÀLCUL D'ERROR - ITERACIÓ SEGÜENT
            !------------------------------------------------
            
            !Estudiem la convergencia en un punt i guardem les dades en un fitxer
            write(unit_file, '(I20, F20.6)') iteracio, T_k_plus_1(punt_i,punt_j)

            !Guardar en el vector conv_p
            conv_p(iteracio) = T_k_plus_1(punt_i,punt_j)

            !Càlcul de l'error
            Error = 0.0d0

            do i = 2, Nx-1
                do j = 2, Ny-1
                    Error = max(Error, abs(T_k_plus_1(i,j) - T_k(i,j)))
                enddo
            enddo

            !Comprovem convergència
            if (Error < precisio) then
                convergencia = .TRUE.
            endif

            !Preparem la següent iteració
            T_k = T_k_plus_1
        
        enddo

        !Guardem el resultat final
        T_final = T_k_plus_1
        
        deallocate(T_k, T_k_plus_1)

        RETURN

    END SUBROUTINE SOLVER_POISSON_2D
