subroutine flow_error
    use global 
    implicit none
    integer :: i
    real :: flow_out, flow_in
    real :: sum
    real :: ai, bi, mi, fa, fb, fm
    
    flow_out = 0.0d0
    
    do i=1,nele
        if(conn(i,5) .eq. -2)then
            ai = xnode(conn(i,2),3);
            bi = xnode(conn(i,4),3);
            mi = xnode(conn(i,3),3);
            fa = nd_vel(conn(i,2),2);
            fb = nd_vel(conn(i,4),2);
            fm = nd_vel(conn(i,3),2);
            
            sum = (bi-ai)/6.0*(fa + 4.0*fm + fb)
            
            flow_out = flow_out + sum
        end if
    end do
    
    flow_in = 0.0d0
    do i=1,nele
        if(conn(i,5) .eq. -1)then
            ai = xnode(conn(i,2),3);
            bi = xnode(conn(i,4),3);
            mi = xnode(conn(i,3),3);
            fa = nd_vel(conn(i,2),2);
            fb = nd_vel(conn(i,4),2);
            fm = nd_vel(conn(i,3),2);
            
            sum = (ai-bi)/6.0*(fa + 4.0*fm + fb)
            
            flow_in = flow_in + sum
        end if
    end do
    
    flow_err = (flow_out-flow_in)/flow_in*100.0
    
end subroutine flow_error