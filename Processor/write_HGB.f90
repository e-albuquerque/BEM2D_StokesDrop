subroutine write_HGB
    use global
    implicit none
    integer :: i,j
    
    open(unit=111, file='mat_H.txt')
    do i=1,2*nnode
        write(111,100)H(i,1:2*nnode)
    end do 
    close(unit=111)
    
    open(unit=112, file='mat_G.txt')
    do i=1,2*nnode
        write(112,100)G(i,1:6*nele)
    end do 
    close(unit=112)
    
    open(unit=113, file='mat_B.txt')
    do i=1,2*nnode
        write(113,100)B(i,nele)
    end do 
    close(unit=113)
    
100 format(200(f8.4))     
end subroutine