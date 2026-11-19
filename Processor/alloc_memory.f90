subroutine alloc_memory
    use global
    implicit none
    
    allocate(xnode(nnode,3),xnode_old(nnode,3),conn(nele,5),bc(nele,13))
    
    allocate(drop_first_last_node(ndrop,3))
    
    allocate(H(2*nnode,2*nnode),G(2*nnode,6*nele),B(2*nnode,nele))
    
    allocate(matA(2*nnode,2*nnode), vector_b(2*nnode), sol(2*nnode), sol_old(2*nnode), d_sol(2*nnode))
    
    allocate(desl_pr(6*nele,6))
    
    allocate(nd_vel(nnode,3), nd_trac(nele,7), nd_vel_old(nnode,3))
end subroutine 