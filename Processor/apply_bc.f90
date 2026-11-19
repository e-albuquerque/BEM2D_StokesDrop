subroutine apply_bc(n_nds, n_el, conect, nds, H, G, B, bc, matA, vector_b)
    implicit none
    integer, intent(in) :: n_nds, n_el
    integer, intent(in) :: conect(n_el,5)
    real, intent(in) :: nds(n_nds,3)
    real, intent(in) :: H(2*n_nds,2*n_nds), G(2*n_nds,6*n_el), B(2*n_nds,n_el)
    real, intent(in) :: bc(n_el,13)
    real, intent(out) :: matA(2*n_nds, 2*n_nds)
    real, intent(out) :: vector_b(2*n_nds)
    
    integer :: el, nd, fd, i, count, de
    integer :: nd3, global_nd, n_des_pr
    integer :: nd_elem(n_nds,5)
    integer :: bc_type
    real :: bc_value
    real :: all_values(6*n_el)
    integer :: share
    integer :: desl_pr(6*n_el,6)
    real :: matB(2*n_nds,6*n_el) 
    integer :: n_no, n_elem, n_fd, n_nd_loc, ind_A, ind_B
    real :: switch(2*n_nds)
    
    nd_elem = 0
    all_values = 0
    matA = 0.0
    count = 0
    desl_pr = 0
    
    do el=1,n_el
        nd3 = conect(el,4)
        do nd=1,3
            global_nd = conect(el,nd+1)
        
            if (nd .eq. 1 .or. nd .eq. 2) then
                nd_elem(global_nd,1) = global_nd
                nd_elem(global_nd,2) = el
                nd_elem(global_nd,3) = nd
                if (nd .eq. 1) then
                    nd_elem(nd3,4)=el
                    nd_elem(nd3,5)=3
                end if
            end if
            do fd=1,2 !freedom degree
                bc_type = nint(bc(el,4*nd-3+2*fd-1))
                bc_value = bc(el,4*nd-2+2*fd-1)
                all_values(6*el+2*nd+fd-8) = bc_value
                if(bc_type .eq. 0)then
                    share = 0
                    i = 1
                    do while (share .eq. 0 .and. i .le. count .and. count .gt. 0)
                        if(global_nd .eq. desl_pr(i,1) .and. fd .eq. desl_pr(i,2))then
                            share = 1
                            desl_pr(i,5) = el
                            desl_pr(i,6) = nd
                        end if
                        i = i + 1
                    end do
                    if(share .eq. 0)then
                        count = count+1
                        desl_pr(count,1)=global_nd
                        desl_pr(count,2)=fd
                        desl_pr(count,3)=el
                        desl_pr(count,4)=nd
                    end if
                end if
            end do
        end do
    end do
    
    n_des_pr = count
    
    matA = H
    matB = G
    
    do de=1,n_des_pr
        n_no = desl_pr(de,1)
        n_fd = desl_pr(de,2)
        n_elem = desl_pr(de,3)
        n_nd_loc = desl_pr(de,4)
        ind_A = (2*n_no-2)+n_fd
        ind_B = (6*n_elem + 2*n_nd_loc - 8)+n_fd
        switch(:) = matB(:,ind_B)
        matB(:,ind_B) = -matA(:,ind_A)
        matA(:,ind_A) = -switch(:)
            if(desl_pr(de,5) .ne. 0)then
                n_elem = desl_pr(de,5)
                n_nd_loc = desl_pr(de,6)
                ind_B = (6*n_elem + 2*n_nd_loc - 8)+n_fd
                matA(:,ind_A)=matA(:,ind_A)-matB(:,ind_B)
                matB(:,ind_B) = 0.0
            end if
    end do
    
    call prod_Axb(2*n_nds, 6*n_el, matB, all_values, vector_b)
    
    do el=1,n_el
        vector_b = vector_b + B(:,el)
    end do
end subroutine 