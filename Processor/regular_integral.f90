subroutine regular_integral(npgr, qsi, w, sf, norm, kappa, radius, r, jacobian, drop, aH, lambda, Ca, h_local, g_local, b_local)
    implicit none
    integer, intent(in) :: npgr, drop
    real, intent(in) :: qsi(npgr), w(npgr), sf(npgr,3), norm(npgr,2), kappa(npgr), radius(npgr), r(npgr,2), jacobian(npgr)
    real, intent(in) :: aH, lambda, Ca
    
    real, intent(out) :: h_local(2,6), g_local(2,6), b_local(2)
    
    real, parameter :: pi = 4.0d0*atan(1.0d0)
    integer :: i
    real :: K11, K12, K22, J11, J12, J22
    real :: rx, ry, r_dot_n, lnr, jw, h_scale, g_scale, b_scale
    
    h_local = 0.0d0
    g_local = 0.0d0
    b_local = 0.0d0
    
    if (drop .eq. 0) then
        h_scale = 1.0d0 / (4.0d0 * pi)
    else
        h_scale = (aH * (1.0d0 - lambda)) / (4.0d0 * pi)
    end if
    g_scale = 1.0d0 / (4.0d0 * pi)
    b_scale = 1.0d0 / (4.0d0 * pi * Ca)
    
    do i = 1, npgr
        rx = r(i,1)
        ry = r(i,2)
        r_dot_n = norm(i,1)*rx + norm(i,2)*ry
        jw = jacobian(i) * w(i)
        
        K11 = (-4.0d0 * rx * rx * r_dot_n / radius(i)) * jw * h_scale
        K12 = (-4.0d0 * rx * ry * r_dot_n / radius(i)) * jw * h_scale
        K22 = (-4.0d0 * ry * ry * r_dot_n / radius(i)) * jw * h_scale
        
        h_local(1,1) = h_local(1,1) + K11 * sf(i,1)
        h_local(1,2) = h_local(1,2) + K12 * sf(i,1)
        h_local(2,1) = h_local(2,1) + K12 * sf(i,1)
        h_local(2,2) = h_local(2,2) + K22 * sf(i,1)
        
        h_local(1,3) = h_local(1,3) + K11 * sf(i,2)
        h_local(1,4) = h_local(1,4) + K12 * sf(i,2)
        h_local(2,3) = h_local(2,3) + K12 * sf(i,2)
        h_local(2,4) = h_local(2,4) + K22 * sf(i,2)
        
        h_local(1,5) = h_local(1,5) + K11 * sf(i,3)
        h_local(1,6) = h_local(1,6) + K12 * sf(i,3)
        h_local(2,5) = h_local(2,5) + K12 * sf(i,3)
        h_local(2,6) = h_local(2,6) + K22 * sf(i,3)
        
        lnr = -log(radius(i))
        J11 = lnr + rx * rx
        J12 = rx * ry
        J22 = lnr + ry * ry
        
        if (drop .eq. 0) then
            J11 = J11 * jw * g_scale
            J12 = J12 * jw * g_scale
            J22 = J22 * jw * g_scale
            
            g_local(1,1) = g_local(1,1) + J11 * sf(i,1)
            g_local(1,2) = g_local(1,2) + J12 * sf(i,1)
            g_local(2,1) = g_local(2,1) + J12 * sf(i,1)
            g_local(2,2) = g_local(2,2) + J22 * sf(i,1)
            
            g_local(1,3) = g_local(1,3) + J11 * sf(i,2)
            g_local(1,4) = g_local(1,4) + J12 * sf(i,2)
            g_local(2,3) = g_local(2,3) + J12 * sf(i,2)
            g_local(2,4) = g_local(2,4) + J22 * sf(i,2)
            
            g_local(1,5) = g_local(1,5) + J11 * sf(i,3)
            g_local(1,6) = g_local(1,6) + J12 * sf(i,3)
            g_local(2,5) = g_local(2,5) + J12 * sf(i,3)
            g_local(2,6) = g_local(2,6) + J22 * sf(i,3)
        else
            b_local(1) = b_local(1) + b_scale * (J11*norm(i,1) + J12*norm(i,2)) * jw * kappa(i)
            b_local(2) = b_local(2) + b_scale * (J12*norm(i,1) + J22*norm(i,2)) * jw * kappa(i)
        end if
        
    end do
    
end subroutine
