subroutine evolution_euler
    use global
    implicit none
    
    integer :: i
    
    xnode_old = xnode
    nd_vel_old = nd_vel
    
    do i=drop_first_last_node(1,2),nnode
        xnode(i,2) = xnode(i,2)+dt*nd_vel(i,2)
        xnode(i,3) = xnode(i,3)+dt*nd_vel(i,3)
    end do
    
    
end subroutine evolution_euler