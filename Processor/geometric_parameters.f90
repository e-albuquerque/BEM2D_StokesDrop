subroutine geometric_parameters(npg, x_nds, sf, dsf, ddsf, xpts, jacobian, norm, tg, curv,  kappa)
    implicit none
    integer, intent(in) :: npg
    real, intent(in) :: x_nds(3,2)
    real, intent(in) :: sf(npg,3), dsf(npg,3), ddsf(npg,3)
    real, intent(in) :: curv(3)
    
    real, intent(out) :: xpts(npg,2), jacobian(npg), norm(npg,2), tg(npg,2), kappa(npg)
    
    real :: curve, a, b, c, radius
    
    real :: ydeltaa, ydeltab, xdeltaa, xdeltab, aslope, bslope, cx, cy, signal
    
    integer :: i
    
    real :: dx(npg,2), ddx(npg,2)
    
    curve = 0.0
    
    do i=1,npg
        xpts(i,:) = sf(i,1)*x_nds(1,:) + sf(i,2)*x_nds(2,:) + sf(i,3)*x_nds(3,:)
    
        dx(i,:) = dsf(i,1)*x_nds(1,:) + dsf(i,2)*x_nds(2,:) + dsf(i,3)*x_nds(3,:)
        
        ddx(i,:) = ddsf(i,1)*x_nds(1,:) + ddsf(i,2)*x_nds(2,:) + ddsf(i,3)*x_nds(3,:)
        
        jacobian(i) = sqrt(dx(i,1)**2 + dx(i,2)**2)
        
        tg(i,:) = dx(i,:)/jacobian(i)
        
        norm(i,1) = tg(i,2)
        norm(i,2) = -tg(i,1)
        
        kappa(i) = sf(i,1)*curv(1) + sf(i,2)*curv(2) + sf(i,3)*curv(3)
    end do
    
end subroutine 