program Prossesing
    use global
    use omp_lib
    implicit none
    integer :: i, j
    real :: t1, t2, t
    !       READING THE PARAMETERS FROM THE PROCESSED MESH (working well)
    
    call read_input 
    
    !finding the first element of the drops
    call first_drop_elem
    
    
    t = 0.0d0
    
    do iterations=1,max_it

        print*, iterations
        
        t1 = omp_get_wtime()
        t = t+dt
        
        do imp_evol = 1,2

            !       CREATING THE MATRICES H, G and B
            print*, 'building_HGB_sing (OpenMP)'
            call building_HGB_sing
            !call write_HGB
            
            !       APPLYING THE BOUNDARY CONDITIONS
            print*,'apply_bc2'
            call apply_bc2 !so far so good
            
            !call write_matA
            
            
            !       SOLVING THE SYSTEM (LAPACK dgesv)
            print*, 'solve_lapack (dgesv)'
            call solve_lapack(2*nnode, matA, vector_b, sol)
            !call Conj_grad(2*nnode, matA, vector_b, sol)
            !call BCGS(2*nnode, matA, vector_b, sol)
            
            !       REORGANIZNG VELOCITIES AND TENSIONS
            print*, 'reorganize'
            call reorganize
            
            !       DROPLET EVOLUTION
            print*,'Droplet evolution'
            if(ndrop .ne. 0)then
                if(imp_evol .eq. 1)then
                    call evolution_euler !first order
                elseif(imp_evol .eq. 2)then
                    call evolution_euler_implicit !second order
                end if
            end if
        end do
        
        !calculate the error associated to the flow rate at the outlet 
        call flow_error
        
        !calc the driver pressure
        call driver_pressure
        
        !writing output data files for post-processing
        call write_output(t)
        
        !       REALOCATING THE NODES
        call nodes_relocation
        call nodes_relocation
        !call nodes_relocation
        !computational time per iteration        
        t2 = omp_get_wtime()
        write(*,*)'time = ', t2 - t1
        
        !conditon to know if the drop is almost in  the end of the channel
        if(xnode(drop_first_last_node(1,2),2) .gt. max_pos)exit
        
    end do
    
    
!100 format(i3,f12.6)
101 format(20(f12.6))    
end program
