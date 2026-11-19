subroutine singular_integral()
    implicit none
    integer, intent(in) :: npgr, node
    real, intent(in) :: qsi(npgr), w(npgr), x_nds(3,2)
    
    real, intent(inout) :: g_local(2,6), b_local(2)
    
    real :: gs(2,6), bs(2)
    real :: gr(2,6), br(2)
    
    real :: pi = 4.0*atan(1.0)
    integer :: i, j
    
    gs = 0.0
    bs = 0.0
    gr = 0.0
    br = 0.0
    
    !Not singular part
    do i=1,npgr
        
    end do
end subroutine