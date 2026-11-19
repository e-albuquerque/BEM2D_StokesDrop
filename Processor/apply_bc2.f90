subroutine apply_bc2
    use global    
    implicit none
    
    integer :: el, nd, fd, i, count, de
    integer :: nd3, global_nd
    integer :: nd_elem(nnode,5)
    integer :: bc_type
    real :: bc_value
    real :: all_values(6*nele)
    integer :: share
    real :: matB(2*nnode,6*nele) 
    integer :: n_no, n_elem, n_fd, n_nd_loc, ind_A, ind_B
    real :: switch(2*nnode)
    
    nd_elem = 0
    all_values = 0.0d0
    matA = 0.0d0
    vector_b = 0.0d0
    count = 0
    desl_pr = 0
    
    do el=1,nele !loop over the elements
        nd3 = conn(el,4) !third node of the element el
        do nd=1,3 !loop over the local nodes
            global_nd = conn(el,nd+1) !global node related to the local node
        
            if (nd .eq. 1 .or. nd .eq. 2) then !case the node nd is the first or second 
                nd_elem(global_nd,1) = global_nd 
                nd_elem(global_nd,2) = el
                nd_elem(global_nd,3) = nd
                if (nd .eq. 1) then !case te node nd is the first one
                    nd_elem(nd3,4)=el !the other element containing that node
                    nd_elem(nd3,5)=3 !it must be the third node of the element
                end if
            end if
            do fd=1,2 !freedom degree
                bc_type = nint(bc(el,4*nd-3+2*fd-1)) !boudary condition type
                bc_value = bc(el,4*nd-2+2*fd-1) !boundary condition value
                all_values(6*el+2*nd+fd-8) = bc_value !saving the boundary condition value on that nd at that fd and in that el
                if(bc_type .eq. 0)then !case it is a velocity boundary condtition
                    share = 0 !turn the share key off
                    i = 1 !counter
                    do while (share .eq. 0 .and. i .le. count .and. count .gt. 0) 
                        if(global_nd .eq. desl_pr(i,1) .and. fd .eq. desl_pr(i,2))then 
                            share = 1 !turn the share key on (this node belongs to other element)
                            desl_pr(i,5) = el !save the el we are in 
                            desl_pr(i,6) = nd !save the local node we are dealing with
                        end if
                        i = i + 1 !counting the iterations
                    end do
                    !basicaly desl_pr looks and saves where the velocity is known! 
                    if(share .eq. 0)then
                        count = count+1 !counting how many predicted velocities do we have
                        desl_pr(count,1)=global_nd !the globe node with predicted veo=locity
                        desl_pr(count,2)=fd !the direction of the known velocity (freedom degree)
                        desl_pr(count,3)=el !the element which this node belongs
                        desl_pr(count,4)=nd !the local node 
                    end if
                end if
            end do
        end do
    end do
    !FOR NOW WE KNOW ALL THE NODES WITH PREDICTED VELOCITIES
    
    n_desl_pr = count !number of predicted velocities
    
    matA = H
    matB = G
    
    do de=1,n_desl_pr !loop over the nodes with predicted velocities
        n_no = desl_pr(de,1) !global node of the predicted velocity node
        n_fd = desl_pr(de,2) !freedom degree which the velocity is known
        n_elem = desl_pr(de,3) !the element containing the node
        n_nd_loc = desl_pr(de,4) !the local number of the node
        ind_A = (2*n_no-2)+n_fd !the indice of the node at the matrix H
        ind_B = (6*n_elem + 2*n_nd_loc - 8)+n_fd !the indice of the node at the matrix G 
        switch(:) = matB(:,ind_B) !saving the column where the velocity is known
        matB(:,ind_B) = -matA(:,ind_A) !passing the column with known velocitie form the right to the left side
        matA(:,ind_A) = -switch(:) !passing the column with unknown velocitie form the left to the right side
            if(desl_pr(de,5) .ne. 0)then !the node belongs to other element 
                n_elem = desl_pr(de,5) !number of the other element the node belongs
                n_nd_loc = desl_pr(de,6) !local number of the node in the other element
                ind_B = (6*n_elem + 2*n_nd_loc - 8)+n_fd !indicates the column of this nd and fd in the matrix G
                matA(:,ind_A)=matA(:,ind_A)-matB(:,ind_B) !pass the unknown traction colunm from the riht to the left 
                matB(:,ind_B) = 0.0 !put zeros on  thet column at the right side
            end if
    end do
    
    call prod_Axb(2*nnode, 6*nele, matB, all_values, vector_b)
    
    !call prod_Axb(2*nnode, 2*nnode, matA, all_values, vector_b)
    
    do el=1,nele
        vector_b = vector_b + B(:,el)
    end do
    
101 format((i3,f9.5))
end subroutine 
