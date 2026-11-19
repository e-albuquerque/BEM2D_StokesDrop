subroutine readmesh
    use global
    implicit none
    
    integer :: i
    
    character (len=30) :: dad !file to be read
    
    dad = trim(inputfilename) // '.txt' !getting the file name taken from "file_name.txt"
    
    open(unit=1,file=dad,status='old') !open the file dad
    
    read(unit=1,fmt=*) !skip a line
    read(unit=1,fmt=*)npts,nseg,ndrop !reading the number of point, n segmnets and n drops (first line)
    
    call alloc_memory1 !alloc memmory for the variables for (xpts, bc_seg, seg, disc, drop_geom)
    
    read(unit=1,fmt=*) !skip a line
    read(unit=1,fmt=*) !skip a line
    read(unit=1,fmt=*)dt,max_pos,max_it !reading the number of point, n segmnets and n drops (first line)
    
    !loop to read the position which define the wall contour
    do i=1,npts
        read(unit=1,fmt=*)xpts(i,1),xpts(i,2),xpts(i,3)
    end do
    
    !loop to read the connection of the segments of the outer contour
    do i=1,nseg
        read(unit=1,fmt=*)seg(i,1),seg(i,2),seg(i,3),seg(i,4),seg(i,5)
    end do
    
    !loop to read the boundary condtions over each segment of the outer contour
    do i=1,nseg
        read(unit=1,fmt=*)bc_seg(i,1),bc_seg(i,2),bc_seg(i,3),bc_seg(i,4),bc_seg(i,5)
    end do
    
    !loop to read the discretization of the segmnets and calculate the total number of elements
    nele=0
    do i=1,nseg
        read(unit=1,fmt=*)disc(i,1),disc(i,2)
        nele = nele+disc(i,2)
    end do
    
    !read the number of elements per drop
    read(unit=1,fmt=*)drop_nele
    !read the viscosity ratio (lambda)
    read(unit=1,fmt=*)lambda
    !read the capillary number
    read(unit=1,fmt=*)Ca
    !read the drops geometries
    do i=1,ndrop
        read(unit=1,fmt=*)drop_geom(i,1), drop_geom(i,2), drop_geom(i,3), drop_geom(i,4)
    end do
    
    !the aspec ratio is always calculate using the diameter of the first drop
    aH = 2.0*drop_geom(1,4)/(xpts(npts,3) - xpts(1,3));
    
    !calculate the total number of elements
    nele = nele + ndrop*drop_nele
    
    !calculating the number of nodes (2 times the number of elements)
    nnode = nele*2
    
    call alloc_memory2 !alloc memmory for the variables for (xnode, bc, conn)
    
    close(unit = 1) !close the open unit
    
100 format(20(f8.5))
end subroutine