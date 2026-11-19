subroutine evolution_euler_implicit
    use global
    implicit none
    
    integer :: i
    
    do i=drop_first_last_node(1,2),nnode
        xnode(i,2) = xnode_old(i,2)+(nd_vel(i,2)+nd_vel_old(i,2))*dt/2.0
        xnode(i,3) = xnode_old(i,3)+(nd_vel(i,3)+nd_vel_old(i,3))*dt/2.0
    end do
    
end subroutine evolution_euler_implicit