subroutine write_matA
    use global
    implicit none
    integer :: i
    
    open(file = 'matA.txt' , unit = 1)
    do i=1,2*nnode
        write(1,100)matA(i,:)
    end do
    close(unit= 1)
    
100 format(100(f8.5))
end subroutine write_matA
