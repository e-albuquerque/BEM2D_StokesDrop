subroutine tecplot_contour
    use global
    implicit none
    
    integer :: i, j
    
    tecplot_file_wall = 'wall_contour.plt'
    tecplot_file_drop = 'drop_contour.plt'
    open(file=tecplot_file_wall , unit = 100)
    open(file=tecplot_file_drop , unit = 105)
    
    if(ndrop .gt. 0)then
        do i=1,drop_first_last_node(1,2)-1
           write(unit=100,fmt=102)xnode(i,2),xnode(i,3)
        end do
        write(unit=100,fmt=102)xnode(1,2),xnode(1,3)
        
        do i=drop_first_last_node(1,2),nnode
            write(unit=105,fmt=102)xnode(i,2),xnode(i,3)
        end do       
        
    else
        do i=1,nnode
           write(unit=100,fmt=102)xnode(i,2),xnode(i,3)
        end do
    end if
     
101 format('zone T=" time=',f8.4,'", i=',i6,', j=',i6,', f=point')
102 format(1x,20(e18.10))  
    
end subroutine 