subroutine prod_Axb(nrows, ncolums, A, b, x)
    implicit none
    integer, intent(in) :: nrows, ncolums
    real, intent(in) :: A(nrows,ncolums), b(ncolums)
    real, intent(out) :: x(nrows)
    
    integer :: i, j
    real :: sum
    
    !$omp parallel do default(shared) private(i, j, sum) schedule(static)
    do i=1,nrows
        sum  = 0.0
        do j=1,ncolums
            sum = sum + A(i,j)*b(j)
        end do
        x(i) = sum
    end do
    !$omp end parallel do
    
end subroutine 