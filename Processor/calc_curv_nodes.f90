subroutine calc_curv_nodes(k, el)
    use global
    implicit none
    real,intent(out) :: k(3)
    integer, intent(in) :: el
    
    integer :: key
    integer :: i
    integer :: signal(3)
    
    real :: x1, y1
    real :: x2, y2
    real :: x3, y3
    real :: x2_b, y2_b !coordinates of the middle node of the element before
    real :: x2_a, y2_a !coordinates of the middle node of the element after
    
    
    k = 0.0d0
    
    x1 = xnode(conn(el,2),2); y1 = xnode(conn(el,2),3);
    x2 = xnode(conn(el,3),2); y2 = xnode(conn(el,3),3);
    x3 = xnode(conn(el,4),2); y3 = xnode(conn(el,4),3);
    !look if the element is the firts or the last of a drop
    key = 0
    do i=1,ndrop
        if(conn(el,2) .eq. drop_first_last_node(i,2))then !the element is the first one of the drop #i
            x2_b = xnode(drop_first_last_node(i,3),2);
            y2_b = xnode(drop_first_last_node(i,3),3);
            x2_a = xnode(conn(el+1,3),2);
            y2_a = xnode(conn(el+1,3),3);
            key = 1
            !write(*,*)conn(el,2),conn(el,3),conn(el,4),drop_first_last_node(i,3)
            exit
        elseif(conn(el,3) .eq. drop_first_last_node(i,3))then !the element is the last one of the drop #i
            x2_b = xnode(conn(el-1,3),2);
            y2_b = xnode(conn(el-1,3),3);
            x2_a = xnode(drop_first_last_node(i,2)+1,2);
            y2_a = xnode(drop_first_last_node(i,2)+1,3);
            key = 1 
            !write(*,*)conn(el,2),conn(el,3),conn(el,4),drop_first_last_node(i,2)+1
            exit
        end if
    end do
    
    if(key .eq. 0)then
        x2_b = xnode(conn(el-1,3),2);
        y2_b = xnode(conn(el-1,3),3);
        x2_a = xnode(conn(el+1,3),2);
        y2_a = xnode(conn(el+1,3),3);
    end if
    
    call calc_circle_curv(k(1), x2_b, y2_b, x1, y1, x2, y2)
    call calc_circle_curv(k(2), x1, y1, x2, y2, x3, y3)
    call calc_circle_curv(k(3), x2, y2, x3, y3, x2_a, y2_a)
    
    call curv_sign(signal(1), x2_b, y2_b, x1, y1, x2, y2)
    call curv_sign(signal(2), x1, y1, x2, y2, x3, y3)
    call curv_sign(signal(3), x2, y2, x3, y3, x2_a, y2_a)
    
    do i=1,3
        k(i) = k(i)*signal(i)
    end do
    
end subroutine calc_curv_nodes