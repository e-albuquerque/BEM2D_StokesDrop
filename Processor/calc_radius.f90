subroutine calc_radius(npg, x_pts, xn, radius, r)
    implicit none
    integer, intent(in) :: npg
    real, intent(in) :: x_pts(npg,2), xn(2)
    
    real, intent(out) :: radius(npg)
    real, intent(out) :: r(npg,2)
    
    integer :: i
    
    do i=1,npg
        r(i,:) = x_pts(i,:) - xn
        
        radius(i) = sqrt(r(i,1)**2 + r(i,2)**2)
         
        r(i,:) = r(i,:)/radius(i)
    end do
end subroutine calc_radius