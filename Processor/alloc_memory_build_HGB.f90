subroutine alloc_memory_build_HGB
    use var_build_HGB
    implicit none
    
    allocate(qsi(npgr),w(npgr))
    
    allocate(eta(npgs),rho(npgs))
    
    allocate(sf(npgr,3),dsf(npgr,3),ddsf(npgr,3))
    
    allocate(sfs_nd1(npgs,3), sfs_nd2p(npgs,3), sfs_nd2n(npgs,3), sfs_nd3(npgs,3))
    allocate(dsfs_nd1(npgs,3), dsfs_nd2p(npgs,3), dsfs_nd2n(npgs,3), dsfs_nd3(npgs,3))
    allocate(ddsfs_nd1(npgs,3), ddsfs_nd2p(npgs,3), ddsfs_nd2n(npgs,3), ddsfs_nd3(npgs,3))
    
    
    allocate(qsis_nd1(npgs), dqsideta_nd1(npgs))
    allocate(qsis_nd2p(npgs), dqsideta_nd2p(npgs))
    allocate(qsis_nd2n(npgs), dqsideta_nd2n(npgs))
    allocate(qsis_nd3(npgs), dqsideta_nd3(npgs))
    
    
    allocate(x_pts(npgr,2), jacobian(npgr), norm(npgr,2), tg(npgr,2), kappa(npgr))
        
    allocate(radius(npgr), r(npgr,2))

    allocate(x_ptss_nd1(npgs,2), jacobians_nd1(npgs), norms_nd1(npgs,2), tgs_nd1(npgs,2), kappas_nd1(npgs))
    allocate(x_ptss_nd2p(npgs,2), jacobians_nd2p(npgs), norms_nd2p(npgs,2), tgs_nd2p(npgs,2), kappas_nd2p(npgs))
    allocate(x_ptss_nd2n(npgs,2), jacobians_nd2n(npgs), norms_nd2n(npgs,2), tgs_nd2n(npgs,2), kappas_nd2n(npgs))
    allocate(x_ptss_nd3(npgs,2), jacobians_nd3(npgs), norms_nd3(npgs,2), tgs_nd3(npgs,2), kappas_nd3(npgs))
    
    allocate(radiuss_nd1(npgs), rs_nd1(npgs,2))
    allocate(radiuss_nd2p(npgs), rs_nd2p(npgs,2))
    allocate(radiuss_nd2n(npgs), rs_nd2n(npgs,2))
    allocate(radiuss_nd3(npgs), rs_nd3(npgs,2))
    
end subroutine 