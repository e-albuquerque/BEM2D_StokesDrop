SUBROUTINE circle_center(xi,yi,xf,yf,r,xc,yc)
    IMPLICIT NONE
    REAL,INTENT(IN) :: xi,yi,xf,yf,r
    REAL,INTENT(OUT) :: xc,yc
    REAL :: xm, ym, b, t1, t2, n1, n2, h
    
    xm = (xi+xf)/2.0; ym = (yi+yf)/2.0;
    b = sqrt((xf-xi)**2.0 + (yf-yi)**2.0);
    t1 = (xf-xi)/b; t2 = (yf-yi)/b;
    n1 = t2; n2 = -t1;
    h = sqrt(abs(r**2.0 - (b/2.0)**2.0));
    
    IF (r .gt. 0.0) THEN
        IF (n1 .eq. 0.0) THEN
            xc = xm;
            yc = ym - n2/abs(n2)*h;
        ELSE    
            xc = -n1/abs(n1)*sqrt(h**2.0*n1**2.0/(n1**2.0 + n2**2.0)) + xm;
            yc = n2/n1*(xc-xm)+ym;
        END IF
    ELSE
        IF (n1 .eq. 0.0) THEN
            xc=xm;
            yc=ym+n2/abs(n2)*h;
        ELSE
            xc=n1/abs(n1)*sqrt(h**2.0*n1**2.0/(n1**2.0+n2**2.0))+xm;
            yc=n2/n1*(xc-xm)+ym;
        END IF
    END IF
    
END SUBROUTINE circle_center