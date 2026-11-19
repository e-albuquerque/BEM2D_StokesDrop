subroutine calc_drop_area(npg, area, x_pts, norm, jacobian, w)
    implicit none
    integer, intent(in) :: npg
    real, intent(in) :: x_pts(npg,2), norm(npg,2), jacobian(npg), w(npg)
    real, intent(inout) :: area
    
    integer :: i
    
    do i=1,npg
        area = area + 0.5*(x_pts(i,1)*norm(i,1) + x_pts(i,2)*norm(i,2))*jacobian(i)*w(i);
    end do
    
end subroutine calc_drop_area