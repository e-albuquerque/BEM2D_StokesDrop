subroutine driver_pressure
    use global
    implicit none
    real :: ai, bi, mi, fa, fb, fm
    integer :: i
    real :: first_y, last_y
    integer :: key
    real :: press_in, press_out
    real :: sum
    
    press_in = 0.0d0
    key = 0
    do i=1,nele
        if(conn(i,5) .eq. -1)then
            ai = xnode(conn(i,2),3);
            if(key .eq. 0)then
                first_y = ai;
                key = 1
            end if
            bi = xnode(conn(i,4),3);
            mi = xnode(conn(i,3),3);
            fa = nd_trac(i,2);
            fb = nd_trac(i,6);
            fm = nd_trac(i,4);
            
            sum = (bi-ai)/6.0*(fa + 4.0*fm + fb)
            
            press_in = press_in + sum
        end if
    end do
    last_y = bi
    press_in = press_in/(last_y-first_y)
    
    
    press_out = 0.0d0
    key = 0
    do i=1,nele
        if(conn(i,5) .eq. -2)then
            ai = xnode(conn(i,2),3);
            if(key .eq. 0)then
                first_y = ai;
                key = 1
            end if
            bi = xnode(conn(i,4),3);
            mi = xnode(conn(i,3),3);
            fa = nd_trac(i,2);
            fb = nd_trac(i,6);
            fm = nd_trac(i,4);
            
            sum = (ai-bi)/6.0*(fa + 4.0*fm + fb)
            
            press_out = press_out + sum
        end if
    end do
    last_y = bi
    press_out = press_out/(last_y-first_y)
    
    d_pressure = press_in-press_out
    
end subroutine driver_pressure