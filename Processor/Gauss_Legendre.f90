subroutine Gauss_Legendre(x1,x2,n,x,w)
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: x1,x2

    real,intent(out) :: x(n),w(n)

    real :: eps,xm,xl,z,p1,p2,p3,pp,z1
    real :: pi = 4.0*atan(1.0);
    integer :: i,j,m

    eps = 3.0**(-14.0);
    m = nint((n+1.0)/2.0);
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);
    !write(*,*)eps,m,xm,xl
    do i=1,m
        z = cos(pi*(i-.25)/(n+.5));
        do while (1 == 1)
            p1=1.0;
            p2=0.0;
            do j=1,n
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1)*p3)/j;
            end do
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
            if (abs(z-z1)<eps) then
                exit
            end if
        end do
        x(i)=xm-xl*z;
        x(n+1-i)=xm+xl*z;
        w(i)=2.0*xl/((1.0-z*z)*pp*pp);
        w(n+1-i)=w(i);
    end do

end subroutine
