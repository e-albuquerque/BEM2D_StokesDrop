module global
    
    character (len=30) :: inputfilename
    character (len=30) :: tecplot_file_wall
    character (len=30) :: tecplot_file_drop
    
    integer :: nnode
    integer :: nele
    integer :: npts
    integer :: nseg
    integer :: ndrop
    integer :: drop_nele
    integer :: drop_first_nd
    
    real :: max_dim
    
    real, allocatable :: xpts(:,:)
    real, allocatable :: bc_seg(:,:)
    real, allocatable :: seg(:,:)
    integer, allocatable :: disc(:,:)
    real,allocatable :: drop_geom(:,:)
    integer,allocatable :: drop_first_last_node(:,:)
    
    
    real :: aH
    real :: lambda
    real :: Ca
    
    real :: dt
    integer :: max_it
    real :: max_pos
    
    real, allocatable :: xnode(:,:)
    real, allocatable :: bc(:,:)
    integer, allocatable :: conn(:,:)
    
end module global 