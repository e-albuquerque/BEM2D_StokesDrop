subroutine apply_bc2
    use global    
    implicit none
    
    integer :: el, nd, fd, i, count, de, j
    integer :: global_nd
    integer :: bc_type
    real :: bc_value
    real :: all_values(6*nele)
    integer :: is_vel_col(6*nele)
    integer :: share
    integer :: n_no, n_elem, n_fd, n_nd_loc, ind_A, ind_B, ind_B2
    
    all_values = 0.0d0
    is_vel_col = 0
    vector_b = 0.0d0
    
    ! Build boundary condition topological map once on the first assembly
    if (iterations .eq. 1 .and. imp_evol .eq. 1) then
        count = 0
        desl_pr = 0
        do el=1,nele
            do nd=1,3
                global_nd = conn(el,nd+1)
                do fd=1,2
                    bc_type = nint(bc(el,4*nd-3+2*fd-1))
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
        n_desl_pr = count
    end if

    do el=1,nele
        do nd=1,3
            do fd=1,2
                all_values(6*el+2*nd+fd-8) = bc(el,4*nd-2+2*fd-1)
            end do
        end do
    end do
    
    matA = H
    
    do de=1,n_desl_pr !loop over nodes with prescribed velocities
        n_no = desl_pr(de,1)
        n_fd = desl_pr(de,2)
        n_elem = desl_pr(de,3)
        n_nd_loc = desl_pr(de,4)
        ind_A = (2*n_no-2)+n_fd
        ind_B = (6*n_elem + 2*n_nd_loc - 8)+n_fd
        is_vel_col(ind_B) = 1

        bc_value = all_values(ind_B)
        if (bc_value .ne. 0.0d0) then
            vector_b(:) = vector_b(:) - H(:, ind_A) * bc_value
        end if

        if(desl_pr(de,5) .ne. 0)then !node is shared between two elements
            n_elem = desl_pr(de,5)
            n_nd_loc = desl_pr(de,6)
            ind_B2 = (6*n_elem + 2*n_nd_loc - 8)+n_fd
            is_vel_col(ind_B2) = 1
            matA(:,ind_A) = -G(:,ind_B) - G(:,ind_B2)
        else
            matA(:,ind_A) = -G(:,ind_B)
        end if
    end do

    ! Accumulate non-zero prescribed tractions (if any)
    do j=1,6*nele
        if (is_vel_col(j) .eq. 0 .and. all_values(j) .ne. 0.0d0) then
            vector_b(:) = vector_b(:) + G(:, j) * all_values(j)
        end if
    end do
    
    ! Accumulate capillary vector B from drop elements
    if (ndrop .gt. 0) then
        do el=drop_f_elem,nele
            vector_b(:) = vector_b(:) + B(:,el)
        end do
    end if
    
101 format((i3,f9.5))
end subroutine
