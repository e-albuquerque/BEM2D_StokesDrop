subroutine sing_reg_integral(x_nds, npgr, qsi, w, sf, norm, kappa, radius, r, jacobian, Ca, deno, qsi0, gns, bns)
    implicit none
    real, intent(in) :: x_nds(3,2)
    integer, intent(in) :: npgr
    real, intent(in) :: qsi(npgr), w(npgr), sf(npgr,3), norm(npgr,2), kappa(npgr)
    real, intent(in) :: radius(npgr), r(npgr,2), jacobian(npgr)
    real, intent(in) :: deno
    integer, intent(in) :: qsi0
    real, intent(in) :: Ca
    
    real, intent(out) :: gns(2,6), bns(2)
    
    integer :: i
    real :: N(2,6), rr(2,2)
    real :: l1, l2, lnrns
    real :: g_st_ns(2,2), b_st_ns(2)
    
    real :: prod(2,6)
    real :: pi
    
    pi = 4.0d0*atan(1.0)
    
    gns = 0.0
    bns = 0.0
    
    do i=1,npgr
        N(1,1) = sf(i,1); N(1,2) = 0.0; N(1,3) = sf(i,2); N(1,4) = 0.0; N(1,5) = sf(i,3); N(1,6) = 0.0;
        N(2,1) = 0.0; N(2,2) = sf(i,1); N(2,3) = 0.0; N(2,4) = sf(i,2); N(2,5) = 0.0; N(2,6) = sf(i,3);
        
        rr(1,1) = r(i,1)*r(i,1); rr(1,2) = r(i,1)*r(i,2);
        rr(2,1) = r(i,2)*r(i,1); rr(2,2) = r(i,2)*r(i,2);
        
        l1 = ((qsi(i)+qsi0)*(x_nds(1,1)-2*x_nds(2,1)+x_nds(3,1)) + x_nds(3,1) - x_nds(1,1))/deno;
        l2 = ((qsi(i)+qsi0)*(x_nds(1,2)-2*x_nds(2,2)+x_nds(3,2)) + x_nds(3,2) - x_nds(1,2))/deno;
        
        lnrns = log(sqrt(l1**2 + l2**2))
        
        g_st_ns(1,1) = -lnrns; g_st_ns(1,2) = 0.0d0
        g_st_ns(2,1) = 0.0d0 ; g_st_ns(2,2) = -lnrns
        
        g_st_ns = (g_st_ns + rr)/(4.0d0*pi)
        
        b_st_ns(1) = g_st_ns(1,1)*norm(i,1) + g_st_ns(1,2)*norm(i,2)
        b_st_ns(2) = g_st_ns(2,1)*norm(i,1) + g_st_ns(2,2)*norm(i,2)
        
        b_st_ns = b_st_ns*kappa(i)/Ca
        
        call mat_prod_2x2_2x6(g_st_ns,N,prod)
        
        gns = gns + prod*jacobian(i)*w(i)
        
        bns = bns + b_st_ns*jacobian(i)*w(i)
        
    end do
    
100 format(10(f7.4))     
end subroutine sing_reg_integral