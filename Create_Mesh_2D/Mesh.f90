!   THIS PROGRAM CREATS A 2D MESH FOR THE A CAPILLARY WITH DROPS INSIDE
!               ------------PARABOLIC SHAPE FUNCTIONS------------
PROGRAM Mesh
    use global    
    IMPLICIT NONE
    
    call getfilename
    call readmesh
    call meshdisc
    call format_bc
    call output_file
    call tecplot_contour
    
    PAUSE
END PROGRAM 