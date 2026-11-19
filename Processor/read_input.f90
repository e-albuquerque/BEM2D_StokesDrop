subroutine read_input
    use global
    implicit none
    integer :: i
    
    inputfilename = 'output.txt'
    
    open(file = inputfilename,  unit = 101) !opennin the file to read the data
    
    !reading the number of nodes, elements and drops
    read(101,*)nnode, nele, ndrop
    
    !alloc memmory for the variables in function os nnode, nele, and ndrop
    call alloc_memory
    
    read(101,*)dt,max_pos,max_it
    
    !skip 2 lines
    read(101,*)
    read(101,*)    
    
    !read the position of the nodes
    do i=1,nnode
        read(101,*)xnode(i,1),xnode(i,2),xnode(i,3)
    end do
    
    !skip 2 lines
    read(101,*)
    read(101,*)
    
    !read the connectivity
    do i=1,nele
        read(101,*)conn(i,1),conn(i,2),conn(i,3),conn(i,4),conn(i,5)
    end do
    
    !skip 2 lines
    read(101,*)
    read(101,*)
    
    !read the boundary conditions
    do i=1,nele
        read(101,*)bc(i,1),bc(i,2),bc(i,3),bc(i,4),bc(i,5),bc(i,6),bc(i,7),bc(i,8),bc(i,9),bc(i,10),bc(i,11),bc(i,12),bc(i,13)
    end do
    
    !skip 2 lines
    read(101,*)
    read(101,*)
    
    if(ndrop .ne. 0)then
        !read the drops firsts and lasts nodes
        do i=1,ndrop
            read(101,*)drop_first_last_node(i,1),drop_first_last_node(i,2),drop_first_last_node(i,3)
        end do
    
        !skip 2 lines
        read(101,*)
        read(101,*)
    
        !read the drop properties    
        read(101,*)lambda, Ca, aH
    end if
    
    close(unit = 101)!closing the file
    
100 format(20(f10.7))
101 format(20(i8))    
    
end subroutine
