SUBROUTINE format_bc
    use global
    IMPLICIT NONE
    
    INTEGER :: s, e, i
    INTEGER :: start_el, end_el, count_el
    INTEGER :: typeBCx1,typeBCy1,typeBCx2,typeBCy2,typeBCx3,typeBCy3
    REAL :: valueBCx1,valueBCy1,valueBCx2,valueBCy2,valueBCx3,valueBCy3
    
    REAL :: mean_vel
    INTEGER :: parab, drop
    REAL :: xc, yc, delta_x, delta_y, p1, p2, xp1, yp1, xp2, yp2
    REAL :: tg(2), norm(2), teta
    
    INTEGER :: nd1, nd2, nd3
    REAL :: xnd1, ynd1, xnd2, ynd2, xnd3, ynd3
    REAL :: vel1, vel2, vel3
    REAL :: d1, d2, d3
    integer :: last_elem
    
    count_el = 0 !element counter
    
    DO s=1,nseg !loop over the segments
        
        start_el = count_el + 1 !start element of the segment
        end_el = count_el + disc(s,2) !end element of the segment
        
        parab = 0 !true false statment to indentify a parabolic velocity profile
        drop = 0 !true false statment to indentify a drop surface condition
        
        
        IF (NINT(bc_seg(s,2)) .eq. 2.0) THEN !the bc is a parabolic velocity profile 
            parab = 1 !bc is a parabolic profile
            IF (seg(s,4) .ne. 0.0) THEN !verifying if the segment is straight
                WRITE(*,*)'Error. The segment must be straight to admit a parabolic velocity profile.'
                PAUSE
                STOP
            ELSE
                typeBCx1 = 0.0;     typeBCx2 = 0.0;     typeBCx3 = 0.0;    
                typeBCy1 = 0.0;     typeBCy2 = 0.0;     typeBCy3 = 0.0;
                mean_vel = bc_seg(s,3); !mean velocity of the parabolic velocity profile
                !verifing the extremes points of the segment
                p1 = seg(s,2);  xp1 = xpts(p1,2);    yp1 = xpts(p1,3);
                p2 = seg(s,3);  xp2 = xpts(p2,2);    yp2 = xpts(p2,3);
                delta_x = xp2 - xp1; delta_y = yp2 - yp1;
                !calculating the tangent and normal unitary vectors of the segment
                !the normal is pointing inside the contour
                tg(1) = delta_x/sqrt(delta_x**2.0 + delta_y**2.0); 
                tg(2) = delta_y/sqrt(delta_x**2.0 + delta_y**2.0);
                norm(1) = -tg(2); norm(2) = tg(1);
                !calculating the center of the segment
                xc = (xp2+xp1)/2.0;  yc = (yp2+yp1)/2.0;        
            END IF
        ELSEIF (NINT(bc_seg(s,2)) .eq. 3.0) THEN !It's the segment of a drop
            drop = 1 !the segment is a drop surface
            IF (seg(s,4) .eq. 0.0) THEN !verifying if the segment is curve
                WRITE(*,*)'Error. The segment must be curved to admit a drop condition on it.'
                PAUSE
                STOP
            ELSE
                typeBCx1 = 1.0;     typeBCx2 = 1.0;     typeBCx3 = 1.0;    
                typeBCy1 = 1.0;     typeBCy2 = 1.0;     typeBCy3 = 1.0;
                valueBCx1 = 0.0;    valueBCx2 = 0.0;    valueBCx3 = 0.0;    
                valueBCy1 = 0.0;    valueBCy2 = 0.0;    valueBCy3 = 0.0;
            END IF
        ELSE !either is a constant velocity condition or a pressure condition 
            typeBCx1 = bc_seg(s,2);     typeBCx2 = bc_seg(s,2);     typeBCx3 = bc_seg(s,2);    
            typeBCy1 = bc_seg(s,4);     typeBCy2 = bc_seg(s,4);     typeBCy3 = bc_seg(s,4);
            valueBCx1 = bc_seg(s,3);    valueBCx2 = bc_seg(s,3);    valueBCx3 = bc_seg(s,3);    
            valueBCy1 = bc_seg(s,5);    valueBCy2 = bc_seg(s,5);    valueBCy3 = bc_seg(s,5);
        END IF        
        
        DO e = start_el, end_el !loop over the elements of the segment
            IF (parab .eq. 1) THEN
                nd1 = conn(e,2); xnd1 = xnode(nd1,2); ynd1 = xnode(nd1,3)
                nd2 = conn(e,3); xnd2 = xnode(nd2,2); ynd2 = xnode(nd2,3)
                nd3 = conn(e,4); xnd3 = xnode(nd3,2); ynd3 = xnode(nd3,3)
                
                d1 = sqrt((xnd1-xc)**2.0 + (ynd1-yc)**2.0)
                d2 = sqrt((xnd2-xc)**2.0 + (ynd2-yc)**2.0)
                d3 = sqrt((xnd3-xc)**2.0 + (ynd3-yc)**2.0)
                
                vel1 = 3.0/2.0*mean_vel*(1.0 - (2.0*d1/sqrt(delta_x**2.0 + delta_y**2.0))**2.0)
                vel2 = 3.0/2.0*mean_vel*(1.0 - (2.0*d2/sqrt(delta_x**2.0 + delta_y**2.0))**2.0)
                vel3 = 3.0/2.0*mean_vel*(1.0 - (2.0*d3/sqrt(delta_x**2.0 + delta_y**2.0))**2.0)
                
                valueBCx1 = vel1*norm(1); valueBCy1 = vel1*norm(2);
                valueBCx2 = vel2*norm(1); valueBCy2 = vel2*norm(2);
                valueBCx3 = vel3*norm(1); valueBCy3 = vel3*norm(2);
            END IF
            
            bc(e,1) = e; bc(e,2) = typeBCx1;  bc(e,3) = valueBCx1;  bc(e,4) = typeBCy1;  bc(e,5) = valueBCy1;
                         bc(e,6) = typeBCx2;  bc(e,7) = valueBCx2;  bc(e,8) = typeBCy2;  bc(e,9) = valueBCy2;
                         bc(e,10) = typeBCx3; bc(e,11) = valueBCx3; bc(e,12) = typeBCy3; bc(e,13) = valueBCy3;
        END DO
        
        count_el = end_el
        
    END DO
    
    last_elem = count_el + 1
    
    do i=1,ndrop
        do e=last_elem,drop_nele+last_elem-1
            bc(e,1) = e; bc(e,2) = 1.0;  bc(e,3) = 0.0;  bc(e,4) = 1.0;  bc(e,5) = 0.0;
                         bc(e,6) = 1.0;  bc(e,7) = 0.0;  bc(e,8) = 1.0;  bc(e,9) = 0.0;
                         bc(e,10) = 1.0; bc(e,11) = 0.0; bc(e,12) = 1.0; bc(e,13) = 0.0;
        end do
        last_elem = e
    end do
            
END SUBROUTINE