subroutine solve_lapack(N, matA, b, x)
    implicit none
    integer, intent(in) :: N
    real(8), intent(inout) :: matA(N,N)
    real(8), intent(in) :: b(N)
    real(8), intent(out) :: x(N)

    integer :: ipiv(N), info
    external :: dgesv

    x = b
    call dgesv(N, 1, matA, N, ipiv, x, N, info)

    if (info .ne. 0) then
        write(*,*) 'Error in LAPACK dgesv: info = ', info
        stop 1
    end if
end subroutine solve_lapack

! OpenMP-parallelized BLAS dgemm override for multi-threaded LAPACK LU factorization (dgetrf)
subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    implicit none
    character*1, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    real(8), intent(in) :: alpha, beta
    real(8), intent(in) :: a(lda,*), b(ldb,*)
    real(8), intent(inout) :: c(ldc,*)
    integer :: i, j, l
    real(8) :: temp

    if (m .le. 0 .or. n .le. 0 .or. ((alpha .eq. 0.0_8 .or. k .le. 0) .and. beta .eq. 1.0_8)) return

    if (alpha .eq. 0.0_8) then
        if (beta .eq. 0.0_8) then
            !$omp parallel do private(i,j) schedule(static) if(m*n > 16384)
            do j = 1, n
                do i = 1, m
                    c(i,j) = 0.0_8
                end do
            end do
        else
            !$omp parallel do private(i,j) schedule(static) if(m*n > 16384)
            do j = 1, n
                do i = 1, m
                    c(i,j) = beta * c(i,j)
                end do
            end do
        end if
        return
    end if

    if (transa .eq. 'N' .or. transa .eq. 'n') then
        if (transb .eq. 'N' .or. transb .eq. 'n') then
            !$omp parallel do private(i,j,l,temp) schedule(static) if(int(m,8)*int(n,8)*int(k,8) > 32768_8)
            do j = 1, n
                if (beta .eq. 0.0_8) then
                    do i = 1, m
                        c(i,j) = 0.0_8
                    end do
                else if (beta .ne. 1.0_8) then
                    do i = 1, m
                        c(i,j) = beta * c(i,j)
                    end do
                end if
                do l = 1, k
                    temp = alpha * b(l,j)
                    if (temp .ne. 0.0_8) then
                        do i = 1, m
                            c(i,j) = c(i,j) + temp * a(i,l)
                        end do
                    end if
                end do
            end do
        else
            !$omp parallel do private(i,j,l,temp) schedule(static) if(int(m,8)*int(n,8)*int(k,8) > 32768_8)
            do j = 1, n
                if (beta .eq. 0.0_8) then
                    do i = 1, m
                        c(i,j) = 0.0_8
                    end do
                else if (beta .ne. 1.0_8) then
                    do i = 1, m
                        c(i,j) = beta * c(i,j)
                    end do
                end if
                do l = 1, k
                    temp = alpha * b(j,l)
                    if (temp .ne. 0.0_8) then
                        do i = 1, m
                            c(i,j) = c(i,j) + temp * a(i,l)
                        end do
                    end if
                end do
            end do
        end if
    else
        if (transb .eq. 'N' .or. transb .eq. 'n') then
            !$omp parallel do private(i,j,l,temp) schedule(static) if(int(m,8)*int(n,8)*int(k,8) > 32768_8)
            do j = 1, n
                do i = 1, m
                    temp = 0.0_8
                    do l = 1, k
                        temp = temp + a(l,i) * b(l,j)
                    end do
                    if (beta .eq. 0.0_8) then
                        c(i,j) = alpha * temp
                    else
                        c(i,j) = alpha * temp + beta * c(i,j)
                    end if
                end do
            end do
        else
            !$omp parallel do private(i,j,l,temp) schedule(static) if(int(m,8)*int(n,8)*int(k,8) > 32768_8)
            do j = 1, n
                do i = 1, m
                    temp = 0.0_8
                    do l = 1, k
                        temp = temp + a(l,i) * b(j,l)
                    end do
                    if (beta .eq. 0.0_8) then
                        c(i,j) = alpha * temp
                    else
                        c(i,j) = alpha * temp + beta * c(i,j)
                    end if
                end do
            end do
        end if
    end if
end subroutine dgemm
