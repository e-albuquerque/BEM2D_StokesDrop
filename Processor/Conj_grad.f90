subroutine Conj_grad(N, matA, b, x)
    implicit none
    integer, intent(in) :: N
    real, intent(in) :: matA(N,N), b(N)
    real, intent(inout) :: x(N)
    real :: r(N), r_st(N), p(N), p_st(N), matAp(N), matAtp_st(N), Ax(N)
    real :: rr_st, matApp_st, alfa, beta, rr_st_n
    integer :: it
    real :: Tol
    
    Tol = 10.0**(-8)
    
    call CG_A_dot_x(N,matA,x,Ax)

    r = b - Ax !resido
    p = r
    r_st = r
    p_st = r_st
    
    it = 1
    
    do while(1 .eq. 1)
        call a_dot_b(N,r,r_st,rr_st)
        call CG_A_dot_x(N,matA,p,matAp)
        call a_dot_b(N,matAp,p_st,matApp_st)
        alfa = rr_st/matApp_st
        x = x + alfa*p
        r = r - alfa*matAp
        matAtp_st = matmul(p_st,matA)
        r_st = r_st - alfa*matAtp_st
        call a_dot_b(N,r,r_st,rr_st_n)
        beta = rr_st_n/rr_st
        p = r + beta*p
        p_st = r_st + beta*p_st
        if(maxval(abs(r)) .lt. Tol)then
            call CG_A_dot_x(N,matA,x,Ax)
            r = b - Ax
             if(maxval(abs(r)) .lt. Tol) exit
        end if
        it = it+1
    end do
    
    !write(*,*)it-1
    
    end subroutine Conj_grad

    
subroutine BCGS(N, matA, b, x) !yet not working well (Wikipedia algorithim )
    implicit none
    integer, intent(in) :: N
    real, intent(in) :: matA(N,N), b(N)
    real, intent(inout) :: x(N)
    real :: r(N), r_st(N), p(N), p_0(N), Ax(N), v_0(N), h(N), s(N), t(N), x_old(N)
    real :: rho, alfa, beta, rho_0, omega_0, ts, tt
    !real :: matAt(N,N)
    integer :: it
    real :: Tol
    
    Tol = 10.0**(-6)
    
    !call tanspose(N,matA,matAt)
    call CG_A_dot_x(N,matA,x,Ax)

    r = b - Ax !resido
    r_st = r
    rho_0 = 1.0d0
    alfa = 1.0d0
    omega_0 = 1.0d0
    v_0  = 0.0d0
    p_0 = 0.0d0
    
    it = 1
    x_old = x
    
    do while(1 .eq. 1)
        call a_dot_b(N,r,r_st,rho)
        beta = (rho/rho_0)/(alfa/omega_0)
        p = r + beta*(p_0 - omega_0*v_0)
        call CG_A_dot_x(N,matA,p,v_0)
        call a_dot_b(N,r_st,v_0,alfa)
        alfa = rho/alfa
        h = x + alfa*p
        if(maxval(abs(alfa*p)) .lt. Tol)then
            x = h
            exit
        end if
        s = r - alfa*v_0
        call CG_A_dot_x(N,matA,s,t)
        call a_dot_b(N,t,s,ts)
        call a_dot_b(N,t,t,tt)
        omega_0 = ts/tt
        x = h + omega_0*s
        if(maxval(abs(x-x_old)) .lt. Tol)exit
        x_old = x
        r = s - omega_0*t
        it = it+1
    end do
    
    write(*,*)it
    
end subroutine BCGS

!       MATRIX VECTOR PRODUCT
subroutine CG_A_dot_x(N,A,x,b)
    implicit none
    integer, intent(in) :: N
    real, intent(in) :: A(N,N), x(N)
    real , intent(out) :: b(N)
    
    b = matmul(A,x)
        
    return
end subroutine CG_A_dot_x

!       VECTORS DOT PRODUCT
subroutine a_dot_b(N,a,b,c)
    implicit none
    integer, intent(in) :: N
    real, intent(in) :: a(N), b(N)
    real, intent(out) :: c
    
    c = dot_product(a,b)
        
    return
end subroutine a_dot_b
    
!       TRANSPOSE MATRIX
subroutine tanspose(N,A,At)
    implicit none
    integer, intent(in) :: N
    real, intent(in) :: A(N,N)
    real, intent(out) :: At(N,N)
    integer :: i,j
    
    do i=1,N
        do j=1,N
            At(i,j) = A(j,i)
        end do
    end do
    
end subroutine
