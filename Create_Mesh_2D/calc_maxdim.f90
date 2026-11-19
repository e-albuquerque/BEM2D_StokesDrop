subroutine calc_maxdim
    use global
    implicit none
    
    real :: max_dl
    
    max_dim = max_dl(npts,nseg,xpts,seg)
    
end subroutine calc_maxdim