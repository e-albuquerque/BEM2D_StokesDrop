subroutine regular_integral(npgr, qsi, w, sf, norm, kappa, radius, r, jacobian, drop, aH, lambda, Ca, h_local, g_local, b_local)
    implicit none
    integer, intent(in) :: npgr, drop
    real, intent(in) :: qsi(npgr), w(npgr), sf(npgr,3), norm(npgr,2), kappa(npgr), radius(npgr), r(npgr,2), jacobian(npgr)
    real, intent(in) :: aH, lambda, Ca
    
    real, intent(out) :: h_local(2,6), g_local(2,6), b_local(2)
    
    real :: pi
    integer :: i
    real :: K_dot_n(2,2), J(2,2)
    real :: r_dot_n
    real :: N(2,6), prod(2,6)
    
    pi = 4.0d0*atan(1.0)
    h_local = 0.0
    g_local = 0.0
    b_local = 0.0
    
    do i = 1, npgr
        r_dot_n = norm(i,1)*r(i,1) + norm(i,2)*r(i,2)
        
        K_dot_n(1,1) = r(i,1)*r(i,1)*r_dot_n
        K_dot_n(1,2) = r(i,1)*r(i,2)*r_dot_n
        K_dot_n(2,1) = r(i,2)*r(i,1)*r_dot_n
        K_dot_n(2,2) = r(i,2)*r(i,2)*r_dot_n
        K_dot_n = -4.0*K_dot_n/radius(i)
        
        J(1,1) = -log(radius(i)) + r(i,1)*r(i,1)
        J(1,2) = r(i,1)*r(i,2)
        J(2,1) = r(i,2)*r(i,1)
        J(2,2) = -log(radius(i)) + r(i,2)*r(i,2)
        
        N(1,1) = sf(i,1); N(1,2) = 0.0; N(1,3) = sf(i,2); N(1,4) = 0.0; N(1,5) = sf(i,3); N(1,6) = 0.0;
        N(2,1) = 0.0; N(2,2) = sf(i,1); N(2,3) = 0.0; N(2,4) = sf(i,2); N(2,5) = 0.0; N(2,6) = sf(i,3);
        
        call mat_prod_2x2_2x6(K_dot_n,N,prod)
        
        if (drop .eq. 0)then
            h_local = h_local + prod*jacobian(i)*w(i)/(4.0*pi)
        else
            h_local = h_local + prod*jacobian(i)*w(i)*(aH*(1.0-lambda))/(4.0*pi)
        end if
        
        call mat_prod_2x2_2x6(J,N,prod)
        
        if (drop .eq. 0) then
            g_local = g_local + prod*jacobian(i)*w(i)/(4.0*pi)
        else
            g_local = 0.0
        end if
        
        if (drop .eq. 0) then
            b_local = 0.0
        else
            b_local(1) = b_local(1) + (1.0/(4.0*pi*Ca))*(J(1,1)*norm(i,1) + J(1,2)*norm(i,2))*jacobian(i)*kappa(i)*w(i)
            b_local(2) = b_local(2) + (1.0/(4.0*pi*Ca))*(J(2,1)*norm(i,1) + J(2,2)*norm(i,2))*jacobian(i)*kappa(i)*w(i)
        end if
        
    end do
    
end subroutine