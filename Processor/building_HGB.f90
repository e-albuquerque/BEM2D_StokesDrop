subroutine building_hgb
    use global
    implicit none
    
    real :: h_local(2,6), g_local(2,6), b_local(2)
    
    integer :: nd, el !counters
    integer :: m,n
    integer :: drop !to identify if the element is on the drop or not
    
    integer :: nd_el(3) !global number of the local nodes
    real :: x_nds(3,2) !position of the local nodes
    real :: xn(2) !position of the fount node
    real, allocatable :: radius(:), r(:,:) !distance of the font node to the int point (both norm and direction)
    real, allocatable :: radiuss(:), rs(:,:) !distance of the font node to the int point (both norm and direction)
    
    !numebr of gaussian points
    integer :: npgr, npgs
    real, allocatable :: qsi(:),w(:)
    real, allocatable :: qsis(:),ws(:)
    
    !shape functions
    real, allocatable :: sf(:,:), dsf(:,:), ddsf(:,:)
    real, allocatable :: sfs(:,:), dsfs(:,:), ddsfs(:,:)
    
    !geometric parameters
    real, allocatable :: x_pts(:,:)
    real, allocatable :: jacobian(:), norm(:,:), tg(:,:), kappa(:)
    real, allocatable :: x_ptss(:,:)
    real, allocatable :: jacobians(:), norms(:,:), tgs(:,:), kappas(:)
    
    !putting zeros on H, G and B
    H = 0.0; G=0.0; B=0.0;
    
    !PARAMTERS FOR THE GAUSSIAN NUMERICAL INTEGRATIONS (BOTH REGULAR AND SINGULAR) (correct)        
    npgr = 6 !it must be even to avoid divisions by zero on the inetgrals
    npgs = 8004 !For now, we are using a large number of gaussian points to calculate the singular integrals (remebring they are weak singulaties)
    
    allocate(qsi(npgr),w(npgr))
    call Gauss_Legendre(-1.0, 1.0, npgr, qsi, w)
    allocate(qsis(npgs),ws(npgs))
    call Gauss_Legendre(-1.0, 1.0, npgs, qsis, ws)
    
    !OBTAINING THE SHAPE FUNCTIONS FOR THE INTEGRATIONS
    allocate(sf(npgr,3),dsf(npgr,3),ddsf(npgr,3))
    call shape_functions(npgr, qsi, sf, dsf, ddsf)
    
    allocate(sfs(npgs,3),dsfs(npgs,3),ddsfs(npgs,3))
    call shape_functions(npgs, qsis, sfs, dsfs, ddsfs)
    
    !ALLOCATING THE GEOMETRIC PARAMETERS
    allocate(x_pts(npgr,2), jacobian(npgr), norm(npgr,2), tg(npgr,2), kappa(npgr))
    allocate(radius(npgr), r(npgr,2))
    
    allocate(x_ptss(npgs,2), jacobians(npgs), norms(npgs,2), tgs(npgs,2), kappas(npgs))
    allocate(radiuss(npgs), rs(npgs,2))
    
    do el=1,nele !loop on the elements
        drop = conn(el,5) !to know if the element is on te drop or not
        !taking the global number of the local nodes (correct)
        nd_el(1) = conn(el,2)
        nd_el(2) = conn(el,3)
        nd_el(3) = conn(el,4)
        !taking the elements nodes position (correct)
        x_nds(1,1) = xnode(nd_el(1),2);
        x_nds(1,2) = xnode(nd_el(1),3);
        x_nds(2,1) = xnode(nd_el(2),2);
        x_nds(2,2) = xnode(nd_el(2),3);
        x_nds(3,1) = xnode(nd_el(3),2);
        x_nds(3,2) = xnode(nd_el(3),3);
        !calculating the geometric parameters over the integration points (correct so far, but without the singular integration points)
        call geometric_parameters(npgr, x_nds, sf, dsf, ddsf, x_pts, jacobian, norm, tg, kappa)
        call geometric_parameters(npgs, x_nds, sfs, dsfs, ddsfs, x_ptss, jacobians, norms, tgs, kappas)
        
        do nd=1,nnode !loop over the nodes
            xn(1) = xnode(nd,2)
            xn(2) = xnode(nd,3)
            !calculating the distance in between the node and the points of the element
            call calc_radius(npgr, x_pts, xn, radius, r)
            call calc_radius(npgs, x_ptss, xn, radiuss, rs)
            
            call regular_integral(npgr, qsi, w, sf, norm, kappa, radius, r, jacobian, drop, aH, lambda, Ca,h_local,g_local,b_local)             
            if (nd .eq. nd_el(1) .or. nd .eq. nd_el(2) .or. nd .eq. nd_el(3)) then !singlar integrals
               call regular_integral(npgs,qsis,ws,sfs,norms,kappas,radiuss,rs,jacobians,drop,aH,lambda,Ca, h_local,g_local,b_local)
            end if
            
            H(nd*2-1:nd*2,nd_el(1)*2-1:nd_el(1)*2) = H(nd*2-1:nd*2,nd_el(1)*2-1:nd_el(1)*2) + h_local(:,1:2)
            H(nd*2-1:nd*2,nd_el(2)*2-1:nd_el(2)*2) = H(nd*2-1:nd*2,nd_el(2)*2-1:nd_el(2)*2) + h_local(:,3:4)
            H(nd*2-1:nd*2,nd_el(3)*2-1:nd_el(3)*2) = H(nd*2-1:nd*2,nd_el(3)*2-1:nd_el(3)*2) + h_local(:,5:6)
        
            G(nd*2-1:2*nd,el*6-5:el*6) = g_local
            
            B(nd*2-1:nd*2,el) = B(nd*2-1:nd*2,el) + b_local
           
        end do !loop over the nodes
        
    end do !loop on the elements
    
    !solving the singulary problem (1/r) *indirect method
    do m=1,nnode
        H(2*m-1:2*m,2*m-1:2*m) = 0.0
        do n=1,nnode
            if (n .ne. m) then
                H(2*m-1:2*m,2*m-1:2*m) = H(2*m-1:2*m,2*m-1:2*m) - H(2*m-1:2*m,2*n-1:2*n)
            end if
        end do
    end do 

100 format(6(f5.3,2x))    
end subroutine 
