subroutine shape_functions_sing
    use var_build_HGB
    implicit none
    
    !real, allocatable :: sfs_nd1(:,:), dsfs_nd1(:,:), ddsfs_nd1(:,:)
    
    !-----FIRST NODE SINGULARITY
    !shape functions for singularity on the first node
    sfs_nd1(1:npgs,1) = 0.5*qsis_nd1(1:npgs)*(qsis_nd1(1:npgs)-1.0)
    sfs_nd1(1:npgs,2) = 1.0-qsis_nd1(1:npgs)**2
    sfs_nd1(1:npgs,3) = 0.5*qsis_nd1(1:npgs)*(qsis_nd1(1:npgs)+1.0)
    
    !first derivate shape functions for singularity on the first node
    dsfs_nd1(1:npgs,1) = qsis_nd1(1:npgs) - 0.5
    dsfs_nd1(1:npgs,2) = -2.0*qsis_nd1(1:npgs)
    dsfs_nd1(1:npgs,3) = qsis_nd1(1:npgs) + 0.5
    
    !second derivate shape functions for singularity on the first node
    ddsfs_nd1(1:npgs,1) = 1.0
    ddsfs_nd1(1:npgs,2) = -2.0
    ddsfs_nd1(1:npgs,3) = 1.0
    
    
    !-----SECOND NODE SINGULARITY
    !shape functions for singularity on the second node, positive segment
    sfs_nd2p(1:npgs,1) = 0.5*qsis_nd2p(1:npgs)*(qsis_nd2p(1:npgs)-1.0)
    sfs_nd2p(1:npgs,2) = 1.0-qsis_nd2p(1:npgs)**2
    sfs_nd2p(1:npgs,3) = 0.5*qsis_nd2p(1:npgs)*(qsis_nd2p(1:npgs)+1.0)
    
    !first derivate shape functions for singularity on the second node, positive segment
    dsfs_nd2p(1:npgs,1) = qsis_nd2p(1:npgs) - 0.5
    dsfs_nd2p(1:npgs,2) = -2.0*qsis_nd2p(1:npgs)
    dsfs_nd2p(1:npgs,3) = qsis_nd2p(1:npgs) + 0.5
    
    !second derivate shape functions for singularity on the second node, positive segment
    ddsfs_nd2p(1:npgs,1) = 1.0
    ddsfs_nd2p(1:npgs,2) = -2.0
    ddsfs_nd2p(1:npgs,3) = 1.0
    
    !shape functions for singularity on the second node, negative segment
    sfs_nd2n(1:npgs,1) = 0.5*qsis_nd2n(1:npgs)*(qsis_nd2n(1:npgs)-1.0)
    sfs_nd2n(1:npgs,2) = 1.0-qsis_nd2n(1:npgs)**2
    sfs_nd2n(1:npgs,3) = 0.5*qsis_nd2n(1:npgs)*(qsis_nd2n(1:npgs)+1.0)
    
    !first derivate shape functions for singularity on the second node, negative segment
    dsfs_nd2n(1:npgs,1) = qsis_nd2n(1:npgs) - 0.5
    dsfs_nd2n(1:npgs,2) = -2.0*qsis_nd2n(1:npgs)
    dsfs_nd2n(1:npgs,3) = qsis_nd2n(1:npgs) + 0.5
    
    !second derivate shape functions for singularity on the second node, negative segment
    ddsfs_nd2n(1:npgs,1) = 1.0
    ddsfs_nd2n(1:npgs,2) = -2.0
    ddsfs_nd2n(1:npgs,3) = 1.0
    
    
    !-----FIRST NODE SINGULARITY
    !shape functions for singularity on the third node
    sfs_nd3(1:npgs,1) = 0.5*qsis_nd3(1:npgs)*(qsis_nd3(1:npgs)-1.0)
    sfs_nd3(1:npgs,2) = 1.0-qsis_nd3(1:npgs)**2
    sfs_nd3(1:npgs,3) = 0.5*qsis_nd3(1:npgs)*(qsis_nd3(1:npgs)+1.0)
    
    !first derivate shape functions for singularity on the third node
    dsfs_nd3(1:npgs,1) = qsis_nd3(1:npgs) - 0.5
    dsfs_nd3(1:npgs,2) = -2.0*qsis_nd3(1:npgs)
    dsfs_nd3(1:npgs,3) = qsis_nd3(1:npgs) + 0.5
    
    !second derivate shape functions for singularity on the third node
    ddsfs_nd3(1:npgs,1) = 1.0
    ddsfs_nd3(1:npgs,2) = -2.0
    ddsfs_nd3(1:npgs,3) = 1.0
    
end subroutine shape_functions_sing