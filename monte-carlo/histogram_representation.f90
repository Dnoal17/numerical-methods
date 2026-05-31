
!-------------------------------
!SUBRUTINA GENERACIÓ HISTOGRAMES
!-------------------------------

!Aquesta subrutina genera un histograma amb ncaixes i calcula els seus errors per representar una densitat de probabilitat
!S'espera que el vector xdata contingui números aleatoris generats a partir d'una densitat de probabilitat g(x) amb algun dels següents mètodes
SUBROUTINE HISTOGRAMA(ndat,xdata,ncaixes,xhisto,histo,errhisto)

    IMPLICIT NONE

    integer :: ncaixes, ndat, i_caixa, i
    double precision :: x_min, x_max, h, p_i
    double precision :: xdata(ndat), xhisto(ncaixes), histo(ncaixes), errhisto(ncaixes)

    !Trobem el mínim i màxim del vector amb tots els números aleatoris
    x_min = minval(xdata)
    x_max = maxval(xdata)

    !Definim una amplada d'interval h a partir d'aquests tenint en compte que volem ncaixes a l'histograma
    h = (x_max-x_min)/ncaixes

    !Inicialitzem tots els outputs
    histo = 0.0d0 
    xhisto = 0.0d0 
    errhisto = 0.0d0

    !Identifiquem a quina caixa pertany cada número aleatori al vector
    do i=1,ndat

        !Si es tracta del cas de xmax hem de tractar-lo a part per no intentar desar-lo a una caixa que no existeix
    	if (xdata(i) .eq. x_max) then
            i_caixa = ncaixes

        !En general
    	else
            i_caixa = int((xdata(i)-x_min)/h) + 1

    	endif

        !Un cop determinada la caixa afegim +1 al contador de punts d'aquella caixa
        histo(i_caixa) = histo(i_caixa) + 1

	enddo

	do i=1,ncaixes

        !Trobem el valor central de cada caixa
        xhisto(i) = x_min + (i-0.5d0)*h

        !Normalitzem l'histograma tenint en compet en número de dades i el tamany de l'interval
        histo(i) = histo(i)/(ndat*h)

        !Calculem l'error a cada caixa tenint en compte que histo(i) segueix una distribució binomial
        p_i = histo(i)*h
        errhisto(i) = dsqrt(p_i*(1.0d0-p_i)/ndat)/h

	enddo

    RETURN

END SUBROUTINE HISTOGRAMA