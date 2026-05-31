
!---------------------
!MÈTODE DE SIMPSOM 3/8
!---------------------

!Aquest és un mètode d'integració numèrica familia dels mètodes de Simpsom
!El nombre d'intèrvals en aplicar el mètode ha de ser MÚLTIPLE DE 3 SÍ O SÍ 
SUBROUTINE METODESIMPSON38(x1, x2, k, funcio, resul)

    IMPLICIT NONE

    double precision :: x1, x2, resul, h, x, suma
    integer :: k, i, N
    double precision, external :: funcio

    !Calculem N com un múltiple de 3
    !Fem 3^k intervals per garantir que sigui múltiple de 3
    !Alternativa: N = 3*2**k seria més flexible
    N = 3 * 2**k
    h = (x2 - x1)/dble(N)
    suma = 0.0d0

    do i = 0, N

        x = x1 + h*i

        !Extrems de l'interval tenen pes 1
        if (i == 0 .or. i == N) then
            suma = suma + funcio(x)

        !Si l'índex és múltiple de 3 i no és extrem, factor 2
        else if (mod(i,3) == 0) then
            suma = suma + 2.0d0 * funcio(x)

        !La resta dels punts tenen factor 3
        else
            suma = suma + 3.0d0 * funcio(x)
        endif

    enddo

    !Multiplicació final pel factor 3h/8
    resul = suma * (3.0d0*h/8.0d0) !Error O(h**4)

    RETURN

END SUBROUTINE METODESIMPSON38
