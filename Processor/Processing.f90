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
    
    t1 = omp_get_wtime()
    write(*,'(A,I0)') 'number of processors = ', omp_get_max_threads()
    
    t = 0.0d0
    
    do iterations=1,max_it

        write(*,'(A,A,I0,A,I0)', advance='no') achar(13), 'time step = ', iterations, '/', max_it
        call flush(6)
        
        t = t+dt
        
        do imp_evol = 1,2

            !       CREATING THE MATRICES H, G and B
            call building_HGB_sing
            
            !       APPLYING THE BOUNDARY CONDITIONS
            call apply_bc2
            
            !       SOLVING THE SYSTEM (LAPACK dgesv)
            call solve_lapack(2*nnode, matA, vector_b, sol)
            
            !       REORGANIZNG VELOCITIES AND TENSIONS
            call reorganize
            
            !       DROPLET EVOLUTION
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
        
        !conditon to know if the drop is almost in  the end of the channel
        if(xnode(drop_first_last_node(1,2),2) .gt. max_pos)exit
        
    end do
    write(*,*)
    t2 = omp_get_wtime()
    write(*,'(A,F0.3,A)') 'processing time = ', t2 - t1, ' s'
    
    
!100 format(i3,f12.6)
101 format(20(f12.6))    
end program
