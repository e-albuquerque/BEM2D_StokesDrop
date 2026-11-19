subroutine alloc_memory1
    use global
    implicit none
    
    allocate(xpts(npts,3)) !i# pt, x, y
    
    allocate(bc_seg(nseg+ndrop,5)) !i# seg, bc type x, bc value x, bc type y, bc value y 
    
    allocate(seg(nseg,5)) !i# seg, 1st node, 2nd node, radius of teh segment, is it inlet or outlet (-1 - in, -2 - out)
    
    allocate(disc(nseg,2)) !i# seg, number of elem
    
    allocate(drop_geom(ndrop,4)) !i# drop, xc, yc, radius
    
    allocate(drop_first_last_node(ndrop,3)) !i# drop, first node, last node
    
end subroutine alloc_memory1