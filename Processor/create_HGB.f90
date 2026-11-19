SUBROUTINE create_HGB(n_nds, n_el, nds, conect, bc, H, G, B)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: n_nds, n_el
    REAL,INTENT(IN) :: nds(n_nds,3)
    REAL,INTENT(IN) :: bc(n_el,13)
    INTEGER,INTENT(IN) :: conect(n_el,4)
    
    REAL,INTENT(OUT) :: H(2*n_nds,2*n_nds)
    REAL,INTENT(OUT) :: G(2*n_nds,6*n_el)
    REAL,INTENT(OUT) :: B(2*n_nds,n_el)
    REAL :: h_local(2,6), g_local(2,6), b_local(2)
    REAL :: h_int(2,6), g_int(2,6), b_int(2)
    
    REAL,ALLOCATABLE :: N(:,:)
    REAL :: mat_N(2,6)
    REAL,ALLOCATABLE :: dN(:,:)
    REAL,ALLOCATABLE :: ddN(:,:)
    
    REAL,ALLOCATABLE :: Ns(:,:)
    REAL,ALLOCATABLE :: dNs(:,:)
    REAL,ALLOCATABLE :: ddNs(:,:)
    
    INTEGER :: npgr, npgs
    REAL,ALLOCATABLE :: qsi(:),w(:)
    REAL,ALLOCATABLE :: eta(:), rho(:)
    
    INTEGER :: i, e, nd, m, nn !counter, element counter, node counter
    
    INTEGER :: nl(3) !local node 
    INTEGER :: n0
    REAL :: x_e(3,2), x0(2)
    !GEOMETRY PARAMETERS - NOT SINGULAR INTEGRATION
    REAL,ALLOCATABLE :: x(:,:)
    REAL,ALLOCATABLE :: tg(:,:), norm(:,:) 
    REAL,ALLOCATABLE :: k(:)
    REAL,ALLOCATABLE :: dxdqsi(:,:), ddxdqsi(:,:)
    REAL,ALLOCATABLE :: dgamadqsi(:), ddgamadqsi(:)
    !GEOMETRY PARAMETERS - SINGULAR INTEGRATION
    REAL,ALLOCATABLE :: x_s(:,:)
    REAL,ALLOCATABLE :: tg_s(:,:), norm_s(:,:) 
    REAL,ALLOCATABLE :: k_s(:)
    REAL,ALLOCATABLE :: dxdqsi_s(:,:), ddxdqsi_s(:,:)
    REAL,ALLOCATABLE :: dgamadqsi_s(:), ddgamadqsi_s(:)
    
    REAL :: r(2), r_mod, r_dl(2), r_dot_n
    
    REAL :: pi = 4.0*atan(1.0)
    
    !       CALCULATING THE POINTS AND WEIGHTS OF THE GAUSS INTEGRATION USING LEGENDRE METHOD
    
    npgr = 8; 
    npgs = 10;
    
    ALLOCATE(qsi(npgr),w(npgr))
    CALL Gauss_Legendre(-1.0, 1.0, npgr, qsi, w) 
    
    ALLOCATE(eta(npgs), rho(npgs))
    CALL Gauss_Logaritimic(npgs, eta, rho) 
    
    !       CALCULATING SHAPE FUNCTIONS AND THE TWO FIRST DERIVATIVES (REGULAR GAUSS INTEGRATION) 
    
    ALLOCATE(N(npgr,3), dN(npgr,3), ddN(npgr,3))
    
    CALL calc_shape_funct(npgr,qsi,N,dN,ddN)
    
    !       CALCULATING SHAPE FUNCTIONS AND THE TWO FIRST DERIVATIVES (LOGARITIMIC GAUSS INTEGRATION) 
    
    ALLOCATE(Ns(npgs,3), dNs(npgs,3), ddNs(npgs,3))
    
    CALL calc_shape_funct(npgs,eta,Ns,dNs,ddNs)
    
    !       ALLOCATING VARIABLES FOR THE GEOMETRY PARAMETERS - Regular integration
    ALLOCATE(x(npgr,2), dxdqsi(npgr,2), ddxdqsi(npgr,2))
    ALLOCATE(dgamadqsi(npgr), ddgamadqsi(npgr))
    ALLOCATE(tg(npgr,2), norm(npgr,2), k(npgr))
    
    !       ALLOCATING VARIABLES FOR THE GEOMETRY PARAMETERS - Singular integration
    ALLOCATE(x_s(npgs,2), dxdqsi_s(npgs,2), ddxdqsi_s(npgs,2))
    ALLOCATE(dgamadqsi_s(npgs), ddgamadqsi_s(npgs))
    ALLOCATE(tg_s(npgs,2), norm_s(npgs,2), k_s(npgs))
    
    !       CALCULATING THE MATRICES H, G and B
    
    H = 0.0;    G = 0.0;    B = 0.0;
    
    DO e=1,n_el !loop over the elements
        
        !coordinates of the element nodes
        nl(1) = conect(e,2);     nl(2) = conect(e,3);     nl(3) = conect(e,4); 
        x_e(1,1) = nds(nl(1),2);    x_e(1,2) = nds(nl(1),3); 
        x_e(2,1) = nds(nl(2),2);    x_e(2,2) = nds(nl(2),3);
        x_e(3,1) = nds(nl(3),2);    x_e(3,2) = nds(nl(3),3);
        
        !calc the geometries
        CALL calc_geom(npgr,N,dN,ddN,x_e,x,dxdqsi,ddxdqsi,dgamadqsi,ddgamadqsi,tg,norm,k)

        !calc the geometries - singular integration
        !CALL calc_geom(npgs,Ns,dNs,ddNs,x_e,x_s,dxdqsi_s,ddxdqsi_s,dgamadqsi_s,ddgamadqsi_s,tg_s,norm_s,k_s)
        
        DO nd = 1,n_nds !loop over the font nodes
            x0(1) = nds(nd,2); x0(2) = nds(nd,3);
            h_local  = 0.0; g_local  = 0.0; b_local  = 0.0;
            
            !NOT SINGULR INTEGRATION
            DO i = 1,npgr !loop over the integration points
                CALL calc_radius(npgr,i,x,x0,norm,r,r_mod,r_dl,r_dot_n)
                
                mat_N(1,1) = N(i,1); mat_N(1,2) = 0.0; mat_N(1,3) = N(i,2); mat_N(1,4) = 0.0; mat_N(1,5) = N(i,3); mat_N(1,6) = 0.0;
                mat_N(2,1) = 0.0; mat_N(2,2) = N(i,1); mat_N(2,3) = 0.0; mat_N(2,4) = N(i,2); mat_N(2,5) = 0.0; mat_N(2,6) = N(i,3);
                
                !SO FAR... THERE IS NO DROPS, JUST RIGID CONTOURS
                
                call fund_solution(r_mod,r_dl,r_dot_n,mat_N,h_int,g_int,b_int)
                
                h_local = h_local + h_int*dgamadqsi(i)*w(i)
                
                g_local = g_local + g_int*dgamadqsi(i)*w(i)               
                
                b_local = b_local + b_int*dgamadqsi(i)*w(i)
            
            END DO
            
            !SINGULAR INTEGARTION
            
            IF (nd .eq. nl(1)) THEN
                
            END IF
            

            H(nd*2-1:nd*2,nl(1)*2-1:nl(1)*2) = H(nd*2-1:nd*2,nl(1)*2-1:nl(1)*2)+h_local(:,1:2)
            H(nd*2-1:nd*2,nl(2)*2-1:nl(2)*2) = H(nd*2-1:nd*2,nl(2)*2-1:nl(2)*2)+h_local(:,3:4)
            H(nd*2-1:nd*2,nl(3)*2-1:nl(3)*2) = H(nd*2-1:nd*2,nl(3)*2-1:nl(3)*2)+h_local(:,5:6)

            G(nd*2-1:nd*2,e*6-5:e*6) = G(nd*2-1:nd*2,e*6-5:e*6) + g_local
            B(nd*2-1:nd*2,e) = B(nd*2-1:nd*2,e) + b_local
                
        END DO
        
        ! Calculo dos termos da diagonal da matriz H (consideraçao de corpo a
        ! temperatura constante)
        DO m = 1,n_nds
            H(2*m-1:2*m,2*m-1:2*m) = 0.0; !zera a diagonal principal
            DO nn = 1, n_nds
                IF(nn .ne. m)THEN
                    H(2*m-1:2*m,2*m-1:2*m) = H(2*m-1:2*m,2*m-1:2*m)- H(2*m-1:2*m,2*nn-1:2*nn);
                END IF
            END DO
        END DO
        
    END DO
    
    OPEN(UNIT = 654, FILE = 'mat_H.txt')
    OPEN(UNIT = 655, FILE = 'mat_G.txt')
    OPEN(UNIT = 656, FILE = 'mat_B.txt')
    
    DO i=1,n_nds*2
        WRITE(654,*)H(i,:)
        WRITE(656,*)B(i,:)
    END DO
    
    DO i=1,2*n_nds
        WRITE(655,*)G(i,:)
    END DO
    
    CLOSE(UNIT = 654)
    CLOSE(UNIT = 655)
    CLOSE(UNIT = 656)
    
    151 FORMAT(f6.2)    
END SUBROUTINE 