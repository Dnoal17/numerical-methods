!-------------------------------------------------------
!MÈTODE DE MONTECARLO CRU (GENERALITZATS PER N-DINENSIÓ)
!-------------------------------------------------------

!Aquesta subrutina realitza integracions de Montecarlo i calcula l'error estimat per una integral de N-dimensions
!El mètode de Montecarlo realitza les integrals aproximant-les al càlcul d'un valor esperat
SUBROUTINE MONTECARLO_CRU(dim,N,fun,a,b,integral,error)

    IMPLICIT NONE

    integer :: k, j, N, dim
    double precision :: a(dim), b(dim), x(N,dim), f(N)
    double precision :: integral, lambda, error, xesperat, moment_2, sigma
    double precision, external :: fun !Funció a integrar

    !Inicialitzem les variables recursives i utilitzem el meu NIUB de llavor
    lambda = 1.0d0
    integral = 0.0d0

    !Fem un bucle per calcular la integral generant N*dim números aleatoris
    do k=1,N

        !Fem un bucle que omple les files de la matriu de DIM nombres aleatoris,
        !cadascun generat uniformement en els seus intervals d'integració corresponents
        do j=1,dim
            x(k,j) = rand()
            x(k,j) = a(j) + (b(j) - a(j))*x(k,j)
        enddo

        !El vector f guarda les imatges de cadascun dels punts
        f(k) = fun(x(k,:))

        !Sumem totes les contribucions
        integral = integral + f(k)

    enddo

    !Determinem a part els factors (b-a) que apareix a la teoria del mètode
    do j=1,dim
        lambda = lambda*(b(j)-a(j))
    enddo

    !Calculem l'aproximació de Montecarlo de la integral
    integral = (lambda/N)*integral

    !Calculem l'esperança i el segon moment
    xesperat = 0.0d0
    moment_2  = 0.0d0

    do k=1,N
        xesperat = xesperat + f(k)
        moment_2  = moment_2  + f(k)**2
    enddo

    xesperat = xesperat / N
    moment_2 = moment_2 / N

    !Desviació típica de la variable aleatòria f(x)
    sigma = dsqrt(moment_2 - xesperat**2)

    !Error del mètode de Montecarlo: σ / √N multiplicat pel volum del domini
    error = (lambda * sigma) / dsqrt(dble(N))

END SUBROUTINE MONTECARLO_CRU
