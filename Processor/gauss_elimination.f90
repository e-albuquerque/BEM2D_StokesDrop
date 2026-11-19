subroutine gauss_elimination(N,mat_A,vet_b,x)
    implicit none
    integer,intent(in) :: N
    real,intent(in) :: mat_A(N,N),vet_b(N)
    real,intent(out) :: x(N)

    integer :: i,j,k,linha,it
    real :: max_c,soma,div
    real :: A(N,N),b(N)
    real :: troca(N),troca_b
    real :: max_d
    integer :: troca_l
    real :: salva_l(N), salva_b

    x = 0.0d0
    A = mat_A;
    b = vet_b;

    !pivotamento
    !do j=1,N
    !    max_d = 0.0
    !    do i=1,N
    !        if (abs(A(i,j)) .gt. max_d) then
    !            max_d = A(i,j)
    !            troca_l = i
    !        end if
    !    end do
    !    salva_l = A(j,:)
    !    salva_b = b(j)
    !    A(j,:) = A(troca_l,:)
    !    b(j) = b(troca_l)
    !    A(troca_l,:) = salva_l
    !    b(troca_l) = salva_b
    !end do

    do k=1,N-1
        do i=k+1,N
            div = (A(i,k)/A(k,k));
            A(i,k:N) = A(i,k:N) - A(k,k:N)*div;
            b(i) = b(i) - b(k)*div;
        end do
    end do

    x(N) = b(N)/A(N,N);
    do i=(N-1),1,-1
        soma = 0.0;
        do j=i+1,N
            soma = soma + A(i,j)*x(j);
        end do
        x(i) = (b(i)-soma)/A(i,i);
    end do

end subroutine gauss_elimination

