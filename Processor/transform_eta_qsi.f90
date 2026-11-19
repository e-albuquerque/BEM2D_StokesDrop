subroutine transform_eta_qsi
    use var_build_HGB
    implicit none
    
    integer :: i
     
    qsis_nd1(1:npgs) = 2.0*eta(1:npgs) - 1.0; dqsideta_nd1(1:npgs) = 2.0
    qsis_nd2p(1:npgs) = eta(1:npgs);          dqsideta_nd2p(1:npgs) = 1.0
    qsis_nd2n(1:npgs) = -eta(1:npgs);         dqsideta_nd2n(1:npgs) = -1.0
    qsis_nd3(1:npgs) = -2.0*eta(1:npgs) + 1.0; dqsideta_nd3(1:npgs) = -2.0
    
end subroutine 