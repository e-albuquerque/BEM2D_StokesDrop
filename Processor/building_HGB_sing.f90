subroutine building_hgb_sing
    use global
    use var_build_HGB
    implicit none
    
    integer :: nd, el !counters
    integer :: m,n
    integer :: drop !to identify if the element is on the drop or not
    
    integer :: nd_el(3) !global number of the local nodes
    real :: x_nds(3,2) !position of the local nodes
    real :: xn(2) !position of the source node
    
    real :: deno
    integer :: qsi0
    real :: curv(3)

    ! Thread-private buffers for OpenMP parallel assembly
    real :: h_loc(2,6), g_loc(2,6), b_loc(2)
    real :: gs_loc(2,6), gns_loc(2,6), bs_loc(2), bns_loc(2)
    real :: x_pts_l(6,2), jacobian_l(6), norm_l(6,2), tg_l(6,2), kappa_l(6)
    real :: radius_l(6), r_l(6,2)
    real :: x_ptss_nd1_l(10,2), jacobians_nd1_l(10), norms_nd1_l(10,2), tgs_nd1_l(10,2), kappas_nd1_l(10)
    real :: x_ptss_nd2p_l(10,2), jacobians_nd2p_l(10), norms_nd2p_l(10,2), tgs_nd2p_l(10,2), kappas_nd2p_l(10)
    real :: x_ptss_nd2n_l(10,2), jacobians_nd2n_l(10), norms_nd2n_l(10,2), tgs_nd2n_l(10,2), kappas_nd2n_l(10)
    real :: x_ptss_nd3_l(10,2), jacobians_nd3_l(10), norms_nd3_l(10,2), tgs_nd3_l(10,2), kappas_nd3_l(10)
    real :: h_el(2*nnode,6)
    
    !putting zeros on H, G and B
    H = 0.0; G=0.0; B=0.0;
    
    !starting drop area calc
    drop_area = 0.0d0
    
    !PARAMTERS FOR THE GAUSSIAN NUMERICAL INTEGRATIONS (BOTH REGULAR AND SINGULAR) (correct)        
    npgr = 6 !it must be even to avoid divisions by zero on the inetgrals
    npgs = 10 !For now, we are using a large number of gaussian points to calculate the singular integrals (remebring they are weak singulaties)
    
    if (iterations .eq. 1 .and. imp_evol .eq. 1) then
        call alloc_memory_build_HGB !alloc all the variables specific for this subroutine 
    end if
    
    call Gauss_Legendre(-1.0, 1.0, npgr, qsi, w)
    
    call Gauss_Logaritimic(npgs, eta, rho)
    
    !OBTAINING THE SHAPE FUNCTIONS FOR THE INTEGRATIONS
    call shape_functions(npgr, qsi, sf, dsf, ddsf)
    
    call transform_eta_qsi !transform the variables from eta to qsi coordinates
    
    call shape_functions_sing

    !$omp parallel do default(shared) &
    !$omp private(el, nd, drop, nd_el, x_nds, xn, deno, qsi0, curv, &
    !$omp         h_loc, g_loc, b_loc, gs_loc, gns_loc, bs_loc, bns_loc, &
    !$omp         x_pts_l, jacobian_l, norm_l, tg_l, kappa_l, radius_l, r_l, &
    !$omp         x_ptss_nd1_l, jacobians_nd1_l, norms_nd1_l, tgs_nd1_l, kappas_nd1_l, &
    !$omp         x_ptss_nd2p_l, jacobians_nd2p_l, norms_nd2p_l, tgs_nd2p_l, kappas_nd2p_l, &
    !$omp         x_ptss_nd2n_l, jacobians_nd2n_l, norms_nd2n_l, tgs_nd2n_l, kappas_nd2n_l, &
    !$omp         x_ptss_nd3_l, jacobians_nd3_l, norms_nd3_l, tgs_nd3_l, kappas_nd3_l, &
    !$omp         h_el) &
    !$omp reduction(+:drop_area) schedule(static)
    do el=1,nele !loop over the elements
        if (conn(el,5) .eq. 1) then
            drop = 1 !to know if the element is on te drop or not
        else 
            drop = 0 !not a drop
        end if
        
        !taking the global number of the local nodes (correct)
        nd_el(1) = conn(el,2)
        nd_el(2) = conn(el,3)
        nd_el(3) = conn(el,4)
        !taking the elements nodes position (correct)
        x_nds(1,1:2) = xnode(nd_el(1),2:3);
        x_nds(2,1:2) = xnode(nd_el(2),2:3);
        x_nds(3,1:2) = xnode(nd_el(3),2:3);
        
        if(drop .eq. 1) then
            call calc_curv_nodes(curv, el)
        else
            curv = 0.0d0
        end if
        
        !calculating the geometric parameters over the integration points
        call geometric_parameters(npgr, x_nds, sf, dsf, ddsf, x_pts_l, jacobian_l, norm_l, tg_l, curv, kappa_l)
        
        call geometric_parameters(npgs, x_nds, sfs_nd1, dsfs_nd1, ddsfs_nd1, x_ptss_nd1_l, jacobians_nd1_l, norms_nd1_l, tgs_nd1_l, curv, kappas_nd1_l)
        call geometric_parameters(npgs, x_nds, sfs_nd2p, dsfs_nd2p, ddsfs_nd2p, x_ptss_nd2p_l, jacobians_nd2p_l, norms_nd2p_l, tgs_nd2p_l, curv, kappas_nd2p_l)
        call geometric_parameters(npgs, x_nds, sfs_nd2n, dsfs_nd2n, ddsfs_nd2n, x_ptss_nd2n_l, jacobians_nd2n_l, norms_nd2n_l, tgs_nd2n_l, curv, kappas_nd2n_l)
        call geometric_parameters(npgs, x_nds, sfs_nd3, dsfs_nd3, ddsfs_nd3, x_ptss_nd3_l, jacobians_nd3_l, norms_nd3_l, tgs_nd3_l, curv, kappas_nd3_l)
        
        do nd=1,nnode !loop over the nodes
            xn(1:2) = xnode(nd,2:3)
            h_loc = 0.0; g_loc = 0.0; b_loc = 0.0;
            
            !calculating the distance in between the node and the points of the element
            call calc_radius(npgr, x_pts_l, xn, radius_l, r_l)
            call regular_integral(npgr, qsi, w, sf, norm_l, kappa_l, radius_l, r_l, jacobian_l, drop, aH, lambda, Ca, h_loc, g_loc, b_loc)
                   
            if (nd .eq. nd_el(1) .or. nd .eq. nd_el(2) .or. nd .eq. nd_el(3)) then !singular integrals
                gs_loc = 0.0d0; gns_loc = 0.0d0
                bs_loc = 0.0d0; bns_loc = 0.0d0
                if(nd .eq. nd_el(1))then
                    deno = 1.0d0
                    qsi0 = -1
                elseif(nd .eq. nd_el(2))then
                    deno = 2.0d0
                    qsi0 = 0
                elseif(nd .eq. nd_el(3))then
                    deno = 1.0d0
                    qsi0 = 1
                end if
                
                call sing_reg_integral(x_nds, npgr, qsi, w, sf, norm_l, kappa_l, radius_l, r_l, jacobian_l, Ca, deno, qsi0, gns_loc, bns_loc)
                call sing_integral(qsi0, Ca, &
                    jacobians_nd1_l, norms_nd1_l, kappas_nd1_l, &
                    jacobians_nd2p_l, norms_nd2p_l, kappas_nd2p_l, &
                    jacobians_nd2n_l, norms_nd2n_l, kappas_nd2n_l, &
                    jacobians_nd3_l, norms_nd3_l, kappas_nd3_l, &
                    gs_loc, bs_loc)
                
                g_loc = gs_loc + gns_loc
                b_loc = bs_loc + bns_loc
                
                if(drop .eq. 1)then
                    g_loc = 0.0
                else
                    b_loc = 0.0
                end if
                
            end if
 
            h_el(nd*2-1:nd*2, 1:6) = h_loc(:, 1:6)
            G(nd*2-1:2*nd, el*6-5:el*6) = g_loc
            B(nd*2-1:nd*2, el) = b_loc
        end do !loop over the nodes

        ! Midpoint node nd_el(2) belongs exclusively to element el
        H(:, nd_el(2)*2-1:nd_el(2)*2) = h_el(:, 3:4)

        ! Endpoint nodes nd_el(1) and nd_el(3) are shared between adjacent elements
        !$omp critical(assemble_H_endpoints)
        H(:, nd_el(1)*2-1:nd_el(1)*2) = H(:, nd_el(1)*2-1:nd_el(1)*2) + h_el(:, 1:2)
        H(:, nd_el(3)*2-1:nd_el(3)*2) = H(:, nd_el(3)*2-1:nd_el(3)*2) + h_el(:, 5:6)
        !$omp end critical(assemble_H_endpoints)
        
        !calc the drop area
        if(drop .eq. 1)then
            call calc_drop_area(npgr, drop_area, x_pts_l, norm_l, jacobian_l, w)
        end if
        
    end do !loop on the elements
    !$omp end parallel do
    
    drop_area = abs(drop_area)
    
    !solving the singulary problem (1/r) *indirect method
    !$omp parallel do default(shared) private(m, n) schedule(static)
    do m=1,nnode
        H(2*m-1:2*m,2*m-1:2*m) = 0.0
        do n=1,nnode
            if (n .ne. m) then
                H(2*m-1:2*m,2*m-1:2*m) = H(2*m-1:2*m,2*m-1:2*m) - H(2*m-1:2*m,2*n-1:2*n)
            end if
        end do
    end do
    !$omp end parallel do 
    
    !calculating the initial area of the drop
    if(iterations .eq. 1 .and. imp_evol .eq. 1)then
        drop_area_ini = drop_area
    end if

100 format(6(f8.5,2x))    
end subroutine 