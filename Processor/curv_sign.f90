subroutine curv_sign(signal, x1, y1, x2, y2, x3, y3)
    implicit none
    integer, intent(out) :: signal
    real, intent(in) :: x1,y1,x2,y2,x3,y3
    
    real :: ydeltaa, xdeltaa, ydeltab, xdeltab
    real :: dxdqsi2, dydqsi2
    real :: n(2)
    real :: aslope, bslope
    real :: cx, cy
    real :: prod
    
    signal = 0
    
    ydeltaa = y2-y1;
    xdeltaa = x2-x1;
    ydeltab = y3-y2;
    xdeltab = x3-x2;
    
    dxdqsi2 = -0.5*x1 + 0.5*x3;
    dydqsi2 = -0.5*y1 + 0.5*y3;
    
    n(1) = 1.0/(sqrt(dxdqsi2**2 + dydqsi2**2))*dydqsi2;
    n(2) = -1.0/(sqrt(dxdqsi2**2 + dydqsi2**2))*dxdqsi2;
    
    aslope = ydeltaa/xdeltaa;
    bslope = ydeltab/xdeltab;
    
    cx = (aslope*bslope*(y1-y3) + bslope*(x1+x2) - aslope*(x2+x3)) / (2.0*(bslope-aslope));
    cy = -1.0 * (cx-(x1+x2)/2.0) / aslope + (y1+y2)/2.0;
    
    prod = (cx-x2)*n(1) + (cy-y2)*n(2);
    if(prod .lt. 0.0)then
        signal = -1
    else 
        signal = 1
    end if
    
end subroutine curv_sign