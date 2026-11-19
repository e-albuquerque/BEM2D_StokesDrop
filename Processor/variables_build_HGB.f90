module var_build_HGB
    
    real :: h_local(2,6), g_local(2,6), b_local(2)
    
    real :: gs(2,6), gns(2,6)
    real :: bs(2), bns(2)
    
    
    !---DISTANCES
    !regular
    real, allocatable :: radius(:), r(:,:) !distance of the font node to the int point (both norm and direction)
    !singular
    real, allocatable :: radiuss_nd1(:), rs_nd1(:,:) !distance of the font node to the int point (both norm and direction)
    real, allocatable :: radiuss_nd2p(:), rs_nd2p(:,:) !distance of the font node to the int point (both norm and direction)
    real, allocatable :: radiuss_nd2n(:), rs_nd2n(:,:) !distance of the font node to the int point (both norm and direction)
    real, allocatable :: radiuss_nd3(:), rs_nd3(:,:) !distance of the font node to the int point (both norm and direction)
    
    !---- GAUSSIAN QUADRATURE -----
    !numebr of gaussian points
    integer :: npgr, npgs
    real, allocatable :: qsi(:),w(:)
    real, allocatable :: eta(:),rho(:)
    
    !singular local coordinates
    real, allocatable :: qsis_nd1(:), dqsideta_nd1(:) !first node of the element
    real, allocatable :: qsis_nd2p(:), dqsideta_nd2p(:) !second node of the element, positive segment
    real, allocatable :: qsis_nd2n(:), dqsideta_nd2n(:) !second node of the element, negative segment
    real, allocatable :: qsis_nd3(:), dqsideta_nd3(:) !third node of the element    
    
    !-------SHAPE FUNCTIONS----------    
    !regular
    real, allocatable :: sf(:,:), dsf(:,:), ddsf(:,:)
    !singulars
    real, allocatable :: sfs_nd1(:,:), dsfs_nd1(:,:), ddsfs_nd1(:,:)
    real, allocatable :: sfs_nd2p(:,:), dsfs_nd2p(:,:), ddsfs_nd2p(:,:)
    real, allocatable :: sfs_nd2n(:,:), dsfs_nd2n(:,:), ddsfs_nd2n(:,:)
    real, allocatable :: sfs_nd3(:,:), dsfs_nd3(:,:), ddsfs_nd3(:,:)
    
    
    !----GEOMETRIC PARAMETERS
    !regular
    real, allocatable :: x_pts(:,:), jacobian(:), norm(:,:), tg(:,:), kappa(:)
    !singular
    real, allocatable :: x_ptss_nd1(:,:), jacobians_nd1(:), norms_nd1(:,:), tgs_nd1(:,:), kappas_nd1(:)
    real, allocatable :: x_ptss_nd2p(:,:), jacobians_nd2p(:), norms_nd2p(:,:), tgs_nd2p(:,:), kappas_nd2p(:)
    real, allocatable :: x_ptss_nd2n(:,:), jacobians_nd2n(:), norms_nd2n(:,:), tgs_nd2n(:,:), kappas_nd2n(:)
    real, allocatable :: x_ptss_nd3(:,:), jacobians_nd3(:), norms_nd3(:,:), tgs_nd3(:,:), kappas_nd3(:)

    !----STATIONARY WALL-TO-WALL INFLUENCE CACHE
    integer :: n_wall_el, n_wall_nd
    real, allocatable :: H_wall_cache(:,:), G_wall_cache(:,:)
    
end module var_build_HGB