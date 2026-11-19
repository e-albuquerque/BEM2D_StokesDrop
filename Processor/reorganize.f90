subroutine reorganize
    use global
    implicit none
    
    real :: bc_value(6*nele)
    real :: switch
    real :: x(2*nnode)
    real :: value_x
    integer :: el, nd_l, fd, nd
    integer :: de
    integer :: n_nd, n_fd, n_el, n_nd_loc
    integer :: ind_x, ind_cd 
    integer :: i
    
    x = sol
    bc_value = 0.0d0
    nd_vel = 0.0d0
    nd_trac = 0.0d0
    
    do el=1,nele
        do nd_l=1,3
            do fd=1,2
                bc_value(6*(el-1)+2*(nd_l-1)+fd) = bc(el,4*nd_l+2*fd-3)
            end do
        end do
    end do
    
    do de=1, n_desl_pr
        n_nd = desl_pr(de,1) 
        n_fd = desl_pr(de,2) 
        n_el = desl_pr(de,3)
        n_nd_loc = desl_pr(de,4)
        ind_x = (2*n_nd - 2) + n_fd
        ind_cd = (6*n_el+2*n_nd_loc-8)+n_fd
        switch = bc_value(ind_cd)
        bc_value(ind_cd) = x(ind_x)
        value_x = x(ind_x)
        x(ind_x) = switch
        if (desl_pr(de,5) .ne. 0) then
            n_el = desl_pr(de,5)
            n_nd_loc = desl_pr(de,6)
            ind_cd = (6*n_el + 2*n_nd_loc - 8) + n_fd
            bc_value(ind_cd) = value_x
        end if
    end do
    
    do nd=1,nnode
        nd_vel(nd,1) = nd
        nd_vel(nd,2) = x(2*(nd-1)+1)
        nd_vel(nd,3) = x(2*nd)
    end do
    
    do el=1,nele
        nd_trac(el,1) = el
        do nd_l=1,3
            do fd=1,2
                nd_trac(el,2*(nd_l-1)+fd+1) = bc_value(6*(el-1)+2*(nd_l-1)+fd)
            end do
        end do
    end do
    
100 format(5(f9.5))
end subroutine reorganize
