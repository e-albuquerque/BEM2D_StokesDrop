REAL FUNCTION max_dl(n_pts,n_seg,pts,seg) !Deteremine the maximum dimension of the problem
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: n_pts, n_seg !numer of points and segments
    REAL,INTENT(IN) :: pts(1:n_pts,1:3) !points defining the contour
    REAL,INTENT(IN) :: seg(1:n_seg,1:5) !segments definition
    
    INTEGER :: i
    REAL :: dl,x1,x2,y1,y2,p1,p2

    max_dl = 0.0;
    DO i = 1,n_seg
       p1 = seg(i,2) ; p2 = seg(i,3);
       x1 = pts(p1,2) ; y1 = pts(p1,3);
       x2 = pts(p2,2) ; y2 = pts(p2,3);
       dl = sqrt((x2-x1)**2.0+(y2-y1)**2.0);
       IF (max_dl .lt. dl) THEN
         max_dl = dl;
       END IF
    END DO
    
END FUNCTION max_dl
