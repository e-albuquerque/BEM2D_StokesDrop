subroutine calc_circle_curv(k, x1, y1, x2, y2, x3, y3)
    implicit none
    real, intent(out) :: k
    real, intent(in) :: x1,y1,x2,y2,x3,y3
    
    real :: a,b,c
    real :: radius
    
    k = 0.0
    
    a = sqrt((x1-x2)**2+(y1-y2)**2);
    b = sqrt((x1-x3)**2+(y1-y3)**2);
    c = sqrt((x3-x2)**2+(y3-y2)**2);
    
    radius = (a*b*c)/sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a));
    
    k = 1.0/radius
end subroutine calc_circle_curv