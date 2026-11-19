SUBROUTINE calc_arc(x1,y1,x2,y2,xc,yc,tet1,tet2)
! Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
! horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
! and the horizontal direction
    IMPLICIT NONE
    REAL,INTENT(IN) :: x1,y1,x2,y2,xc,yc
    REAL,INTENT(OUT) :: tet1, tet2
    REAL :: dx1,dx2,dy1,dy2,pi

    pi = 4.0d0*atan(1.0d0);

    dx1 = x1 - xc; dy1 = y1 - yc;
    dx2 = x2 - xc; dy2 = y2 - yc;

    ! Computation of tet1
    IF (dy1 .eq. 0.0) THEN				!% The point 1 and the center have the same y coordinate
        IF (x1 .gt. xc) THEN
          tet1 = 0.0;
        ELSE  !% (x1 < xc)
          tet1 = pi;
        END IF
      ELSE IF (dx1 .eq. 0.0) THEN !				% The point 1 and the center have the same x coordinate
        IF (y1 .gt. yc) THEN
          tet1 = pi/2.0;
        ELSE  ! (y1 < yc)
          tet1 = -pi/2.0;
        END IF
      ELSE  !(dx1~=0 e dy1~=0)
        tet1 = atan(dy1/dx1);
        IF (dx1 .lt. 0.0 .and. tet1 .lt. 0.0) THEN
          tet1 = pi + tet1;
        ELSE IF (dx1 .lt. 0.0 .and. tet1 .gt. 0.0) THEN
          tet1 = -pi + tet1;
        END IF
      END IF

! Computation of tet2
      IF (dy2 .eq. 0.0) THEN				! The point 2 and the center have the same y coordinate
        IF (x2 .gt. xc) THEN
          tet2 = 0.0;
        ELSE  ! (x2 < xc)
          tet2 = pi;
        END IF
      ELSE IF (dx2 .eq. 0.0)  THEN			!% The point 2 and the center have the same x coordinate
        IF (y2 .gt. yc) THEN
          tet2 = pi/2.0;
        ELSE  !% (y2 < yc)
          tet2 = -pi/2.0;
        END IF
      ELSE  ! (dx2/=0 e dy2/=0)
        tet2 = atan(dy2/dx2);
        IF (dx2 .lt. 0.0 .and. tet2 .lt. 0.0) THEN
          tet2 = pi + tet2;
        ELSE IF (dx2 .lt. 0.0 .and. tet2 .gt. 0.0) THEN
          tet2 = -pi + tet2;
        END IF
      END IF

END SUBROUTINE calc_arc
