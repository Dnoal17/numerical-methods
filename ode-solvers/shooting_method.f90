!-----------------------------------
!MÈTODE DE TIR + MÈTODE DE LA SECANT
!-----------------------------------

!A continuació es presenta un exemple d'implementació dels mètodes de tir amb la secant per intgrar l'equació d'Schrodinger independent del temps
!Aquestes subrutines no m'ha estat possible generalitzar-les per qualsevol funció
!Modificar-la per integrar altres casos

subroutine metodo_tiro_secante(E1,E2,n_passos,phi,eps)
    
    implicit none

    integer :: n_passos, i
    double precision :: E1, E2, E3, phi_L_1, phi_L_2, phi_L_3
    double precision :: phi(n_passos), dx, eps 

    double precision :: E, Pi, a, b, L, beta, phi0, dphi0
    COMMON /CONDICIONS_INICIALS/ phi0, dphi0
    COMMON /PARAMETRES/ E, Pi, a, b, L, beta

    !Calculem el pas
    dx = (b - a) / (n_passos - 1)

    !Evaluem E en E1 i trobem phi(L)
    !Per estalviar linies de codi utitlizem una subrutina per fer la integració
    E = E1
    call integrar_schrodinger(dx, n_passos, phi_L_1, phi)

    !Evaluem E en E2 i trobem phi(L)
    !Per estalviar linies de codi utitlizem una subrutina per fer la integració
    E = E2
    call integrar_schrodinger(dx, n_passos, phi_L_2, phi)

    !Fixem 100 iteracions màximes per convergir per evitar problemes i utitlizem el mètode de la secant
    !Escriurem al document dins la subrutina per poder saber com convergeix el mètode de la secant
    write(10,*) "ITERACIÓ SECANT, VALOR D'ENERGIA"
    
    do i = 1, 100

        !Si els valors són molt pròxims parem el bucle (divisió entre 0)
        if (abs(phi_L_2 - phi_L_1) < 1.d-14) then
            exit
        endif

        !Calculem la nova energia E3 en la intersecció de la secant amb les abcisses
        E3 = (E1*phi_L_2 - E2*phi_L_1)/(phi_L_2 - phi_L_1)

        !Evaluem E en E3 i trobem phi(L)
        E = E3
        call integrar_schrodinger(dx, n_passos, phi_L_3, phi)

        !Escrivim al fitxer la convergència del valor de E
        write(10,'(I10,F16.9)') i, E

        !Si convergeix segons la tolerància acceptada (eps) sortim acabem el mètode i la variable queda guardada al COMMON (E=E3)
        if (abs(phi_L_3) < eps) then

            write(10,'(/)')

            exit

        endif

        !Si en aquest pas no ha convergit inicialitzem els nous valors i tornem a començar
        E1 = E2
        E2 = E3
        phi_L_1 = phi_L_2
        phi_L_2 = phi_L_3

    end do

    return

end subroutine metodo_tiro_secante

!--------------------------------------------------------
!INTEGRACIÓ L'EQUACIÓ D'SCHRODINGER DINS EL MÈTODE DE TIR
!--------------------------------------------------------

subroutine integrar_schrodinger(dx, n_passos, phi_L, phi_x)
    
    implicit none
    
    integer :: i, n_passos
    double precision :: dx, x, phi_L
    double precision :: y(2), y_actualitzada(2), phi_x(n_passos)

    double precision :: E, Pi, a, b, L, beta, phi0, dphi0
    COMMON /CONDICIONS_INICIALS/ phi0, dphi0
    COMMON /PARAMETRES/ E, Pi, a, b, L, beta

    !Inicialitzem els vectors amb les condicions inicials
    y(1)=phi0
    y(2)=dphi0
    phi_x(1)=y(1)
    x = a
    !Integrem l'ecuació de Schrodinger independent del temps utilitzant la subrutina del mètode Rk4
    !Aquesta subrutina està programada per calcular la funció d'ona només al següent pas temporal
    !Utitlizem un vector phi_x(n_passos) que recull el valor de la funció d'ona en cada pas
    do i=2, n_passos

        !Ens ajudem d'un vector auxiliar y_actualitzada per poder utitlitzar RK4 adientment
        call miRungeKutta4(x, dx, y, 2, y_actualitzada)
        y = y_actualitzada
        phi_x(i)=y(1)

        !Actualitzem la posicio
        x = a + dx*(i-1)

    end do

    !Phi(L) està guardat a l'última entrada del evctor phi_x
    phi_L = phi_x(n_passos)

end subroutine integrar_schrodinger
