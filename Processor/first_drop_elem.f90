subroutine first_drop_elem
    use global
    implicit none
    
    integer :: el
    
    do el=1,nele
        if(conn(el,5) .eq. 1)then
            drop_f_elem = el !the first drop elem
            exit
        end if
    end do
    
end subroutine first_drop_elem