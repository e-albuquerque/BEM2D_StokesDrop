subroutine sing_integral(qsi0, Ca, &
    jac_nd1, nrm_nd1, kap_nd1, &
    jac_nd2p, nrm_nd2p, kap_nd2p, &
    jac_nd2n, nrm_nd2n, kap_nd2n, &
    jac_nd3, nrm_nd3, kap_nd3, &
    gs_out, bs_out)
    use var_build_HGB
    implicit none
    integer, intent(in) :: qsi0
    real, intent(in) :: Ca
    real, intent(in) :: jac_nd1(npgs), nrm_nd1(npgs,2), kap_nd1(npgs)
    real, intent(in) :: jac_nd2p(npgs), nrm_nd2p(npgs,2), kap_nd2p(npgs)
    real, intent(in) :: jac_nd2n(npgs), nrm_nd2n(npgs,2), kap_nd2n(npgs)
    real, intent(in) :: jac_nd3(npgs), nrm_nd3(npgs,2), kap_nd3(npgs)
    real, intent(out) :: gs_out(2,6), bs_out(2)
    integer :: i
    real :: N1(2,6), N2(2,6)
    real :: pi
    
    pi = 4.0d0*atan(1.0d0)
    
    gs_out = 0.0d0
    bs_out = 0.0d0
    
    do i = 1, npgs
        if(qsi0 .eq. -1)then
             N1(1,1) = sfs_nd1(i,1); N1(1,2) = 0.0; N1(1,3) = sfs_nd1(i,2); N1(1,4) = 0.0; N1(1,5) = sfs_nd1(i,3); N1(1,6) = 0.0;
             N1(2,1) = 0.0; N1(2,2) = sfs_nd1(i,1); N1(2,3) = 0.0; N1(2,4) = sfs_nd1(i,2); N1(2,5) = 0.0; N1(2,6) = sfs_nd1(i,3);    
             
             gs_out = gs_out + N1*jac_nd1(i)*rho(i)*dqsideta_nd1(i)/(4.0*pi)
             
             bs_out = bs_out + 1.0/(4.0*pi*Ca)*kap_nd1(i)*nrm_nd1(i,:)*jac_nd1(i)*rho(i)*dqsideta_nd1(i)
            
        elseif(qsi0 .eq. 0)then
            N1(1,1) = sfs_nd2p(i,1); N1(1,2) = 0.0; N1(1,3) = sfs_nd2p(i,2); N1(1,4) = 0.0; N1(1,5) = sfs_nd2p(i,3); N1(1,6) = 0.0;
            N1(2,1) = 0.0; N1(2,2) = sfs_nd2p(i,1); N1(2,3) = 0.0; N1(2,4) = sfs_nd2p(i,2); N1(2,5) = 0.0; N1(2,6) = sfs_nd2p(i,3);
            
            N2(1,1) = sfs_nd2n(i,1); N2(1,2) = 0.0; N2(1,3) = sfs_nd2n(i,2); N2(1,4) = 0.0; N2(1,5) = sfs_nd2n(i,3); N2(1,6) = 0.0;
            N2(2,1) = 0.0; N2(2,2) = sfs_nd2n(i,1); N2(2,3) = 0.0; N2(2,4) = sfs_nd2n(i,2); N2(2,5) = 0.0; N2(2,6) = sfs_nd2n(i,3);
            
            gs_out = gs_out + N1*jac_nd2p(i)*rho(i)*dqsideta_nd2p(i)/(4.0*pi) - N2*jac_nd2n(i)*rho(i)*dqsideta_nd2n(i)/(4.0*pi)
            
            bs_out = bs_out + 1.0/(4.0*pi*Ca)*kap_nd2p(i)*nrm_nd2p(i,:)*jac_nd2p(i)*rho(i)*dqsideta_nd2p(i) - &
&                             1.0/(4.0*pi*Ca)*kap_nd2n(i)*nrm_nd2n(i,:)*jac_nd2n(i)*rho(i)*dqsideta_nd2n(i)
            
        elseif(qsi0 .eq. 1)then
            N1(1,1) = sfs_nd3(i,1); N1(1,2) = 0.0; N1(1,3) = sfs_nd3(i,2); N1(1,4) = 0.0; N1(1,5) = sfs_nd3(i,3); N1(1,6) = 0.0;
            N1(2,1) = 0.0; N1(2,2) = sfs_nd3(i,1); N1(2,3) = 0.0; N1(2,4) = sfs_nd3(i,2); N1(2,5) = 0.0; N1(2,6) = sfs_nd3(i,3);
            
            gs_out = gs_out - N1*jac_nd3(i)*rho(i)*dqsideta_nd3(i)/(4.0*pi)
            
            bs_out = bs_out - 1.0/(4.0*pi*Ca)*kap_nd3(i)*nrm_nd3(i,:)*jac_nd3(i)*rho(i)*dqsideta_nd3(i)
        end if
        
    end do
    
101 format(10(f8.5))    
end subroutine sing_integral 
