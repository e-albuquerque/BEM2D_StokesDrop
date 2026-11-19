subroutine shape_functions(n,x,sf,dsf,ddsf)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: x(n)
    real, intent(out) :: sf(n,3), dsf(n,3), ddsf(n,3)
    
    integer :: i
    
    do i=1,n
        sf(i,1) = x(i) * (x(i) - 1.0) / 2.0
        sf(i,2) = (1.0 - x(i)) * (1.0 + x(i))
        sf(i,3) = x(i) * (x(i) + 1.0) / 2.0
        
        dsf(i,1) = x(i) - 0.5
        dsf(i,2) = - 2.0 * x(i)
        dsf(i,3) = x(i) + 0.5
        
        ddsf(i,1) = 1.0
        ddsf(i,2) = -2.0
        ddsf(i,3) = 1.0
    end do
    
end subroutine