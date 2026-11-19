module global
    
    character (len=30) :: inputfilename
    character (len=30) :: tecplot_file_wall
    character (len=30) :: tecplot_file_drop
    character (len=30) :: tecplot_file_press
    character (len=30) :: tecplot_file_flowerr
    character (len=30) :: tecplot_file_areaerr
      
    integer :: nnode
    integer :: nele
    integer :: ndrop
    integer :: npts
    integer :: nseg
    
    
    integer :: iterations !take care to not use it in any subroutine
    integer :: imp_evol
    
    real :: max_dim
    real :: flow_err, d_pressure
    real :: drop_area, drop_area_ini
    
    real, allocatable :: xpts(:,:)
    real, allocatable :: bc_seg(:,:)
    real, allocatable :: seg(:,:)
    integer, allocatable :: disc(:,:)
    
    real :: aH
    real :: lambda
    real :: Ca
    real :: dt
    real :: max_pos
    integer :: max_it
    
    real, allocatable :: xnode(:,:)
    real, allocatable :: xnode_old(:,:)
    real, allocatable :: bc(:,:)
    integer, allocatable :: conn(:,:)
    
    real, allocatable :: H(:,:), G(:,:), B(:,:)
    real, allocatable :: matA(:,:), vector_b(:), sol(:), sol_old(:), d_sol(:)
    
    integer, allocatable :: desl_pr(:,:)
    integer :: n_desl_pr
    
    real,allocatable :: nd_vel(:,:), nd_trac(:,:)
    real,allocatable :: nd_vel_old(:,:)
    
    integer,allocatable :: drop_first_last_node(:,:)
    integer :: drop_f_elem
    
end module global