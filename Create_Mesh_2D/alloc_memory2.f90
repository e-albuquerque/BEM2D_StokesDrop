subroutine alloc_memory2
    use global
    implicit none
    
    allocate(xnode(nnode,3)) !i# nd, x, y
    
    allocate(bc(nele,13)) !bc for each node of the element
    
    allocate(conn(nele,5)) !i# elem, 1st nd, 2nd nd, 3th nd, type of element (0 = wall, 1 = drop, -1 = inlet, -2 = outlet)
    
end subroutine alloc_memory2