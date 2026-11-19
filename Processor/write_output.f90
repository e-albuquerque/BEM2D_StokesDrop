subroutine write_output(t)
    use global
    implicit none
    real, intent(in) :: t
    
    integer :: i
    
    out_file_wall = 'wall_contour.dat'
    out_file_press = 'pressure.dat'
    out_file_flowerr = 'flow_error.dat'
    out_file_areaerr = 'area_error.dat'
    
    if (iterations .eq. 1) then
        open(file=out_file_press, unit=110, status='replace')
        open(file=out_file_flowerr, unit=115, status='replace')
        open(file=out_file_areaerr, unit=120, status='replace')
    else
        open(file=out_file_press, unit=110, position='append')
        open(file=out_file_flowerr, unit=115, position='append')
        open(file=out_file_areaerr, unit=120, position='append')
    end if
    
    if (iterations .eq. 1 .or. mod(iterations, 10) .eq. 0) then
        write(out_file_drop,"(A,i5,A)") "drop_t=", nint(t*1000), ".dat"
        open(file=out_file_drop, unit=105, status='replace')
        open(file=out_file_wall, unit=100, status='replace')
        if (ndrop .gt. 0) then
            do i=1,drop_first_last_node(1,2)-1
                write(unit=100,fmt=102) xnode(i,2), xnode(i,3)
            end do
            write(unit=100,fmt=102) xnode(1,2), xnode(1,3)

            do i=drop_first_last_node(1,2),nnode
                write(unit=105,fmt=102) xnode(i,2), xnode(i,3)
            end do
            write(unit=105,fmt=102) xnode(drop_first_last_node(ndrop,2),2), xnode(drop_first_last_node(ndrop,2),3)
        else
            do i=1,nnode
                write(unit=100,fmt=102) xnode(i,2), xnode(i,3)
            end do
        end if
        close(unit=100)
        close(unit=105)
    end if
    
    write(unit=110, fmt=102) t, d_pressure
    write(unit=115, fmt=102) t, flow_err
    write(unit=120, fmt=102) t, (drop_area-drop_area_ini)/(drop_area_ini)*100.0
    
    close(unit=110)
    close(unit=115)
    close(unit=120)
    
102 format(1x,20(e18.10))
500 format(2(f11.6))
    
end subroutine write_output
