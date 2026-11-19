SUBROUTINE create_nodes_conn
    use global
    IMPLICIT NONE
    
    INTEGER :: n,s,e,ce,i,j,k !node counter, segment counter, local(segment) element, global element counter, auxiliar counter 
    INTEGER :: p1,p2 !start and end point of the segment
    REAL :: xp1, yp1, xp2, yp2 !coordinates of the start and end points of the segment
    INTEGER :: nd1, nd2, nd3 !local nodes of the elements to creat ethe conectivity matrix
    REAL :: x1, y1, x2, y2, x3, y3 !coordinates of the local nodes
    
    REAL :: delta_x, delta_y !segment size
    REAL :: dx, dy !distance between nodes
    REAL :: r, xc, yc, r1, r2 !radius, center coordinates of the arc
    REAL :: tet1, tet2, tet, dtet !angles related to the arc
    INTEGER :: sig !location of the center in relation to the segment sig = 1 -> left sig = -1 -> right
    
    INTEGER :: ini_pt, first_pt, first_nd !key to know if the point is starting a new contour, point and node number of this point
    
    REAL :: pi = 4.0d0*atan(1.0d0) !pi 3.14...
    
    ! CALCULATING THE MESH NODES COORDINATES
    
    ini_pt = 1 !true or false statement to know if a new contour will start 
    n = 1 !noudes counter
    ce = 1 !elements counter
    
    DO s = 1,nseg !loop over the segments on the walls
        
        p1 = NINT(seg(s,2)); xp1 = xpts(p1,2); yp1 = xpts(p1,3) !first segment point
        p2 = NINT(seg(s,3)); xp2 = xpts(p2,2); yp2 = xpts(p2,3) !second segment point
        
        IF (abs(seg(s,4)) .lt. 0.001) THEN !the segment is a straight line
            
            delta_x = xp2-xp1; delta_y = yp2-yp1; !segment length
            
            dx = delta_x/(2.0*disc(s,2)); dy = delta_y/(2.0*disc(s,2)); !distance between nodes
                    
        ELSE !the segment is an arc
            
            r = seg(s,4) !segment radius
            
            CALL circle_center(xp1,yp1,xp2,yp2,r,xc,yc) !arc center
            
            !verifying the input data
            r1 = sqrt((xc-xp1)**2 + (yc-yp1)**2) 
            r2 = sqrt((xc-xp2)**2 + (yc-yp2)**2)
            
            IF (abs(r1-r2) .lt. max_dim*0.00001) THEN
            
                CALL calc_arc(xp1,yp1,xp2,yp2,xc,yc,tet1,tet2) !points angles in relation to the horizontal
                
                !angle in between the segments of the arc
                IF (tet2 .lt. tet1) THEN
                    tet2 = tet2 + 2.0*pi
                END IF
                
                IF (seg(s,4) .gt. 0.0) THEN
                    tet = abs(tet2 - tet1)
                    sig = 1
                ELSE
                    tet = 2.0*pi - abs(tet2 - tet1)
                    sig = -1
                END IF
                
                !angle distance between the nodes
                dtet = tet/(2.0*disc(s,2))
                
            ELSE
                
                WRITE(*,*)'Error! The contour is not valid!'
                PAUSE
                EXIT
                
            END IF
        END IF
        
        DO e=1,disc(s,2)
            
            !INICIATING THE CONTOUR
            IF (ini_pt .eq. 1) THEN
                
                xnode(n,1) = n; 
                xnode(n,2) = xp1; 
                xnode(n,3) = yp1;
                nd1 = n
                first_nd = nd1
                n = n + 1
                ini_pt = 0
                first_pt = p1        
                x1 = xp1; y1 = yp1;
                
            ELSE
                
                nd1 = n - 1   
                
            END IF
            
            IF (seg(s,4) .eq. 0.0) THEN !segment is a straight segment
              
                x2 = x1 + dx;   y2 = y1 + dy;
                xnode(n,1) = n; xnode(n,2) = x2; xnode(n,3) = y2;
                nd2 = n
                n = n + 1
                
                IF (p2 .eq. first_pt .and. e .eq. disc(s,2)) THEN !this elemnt is the last one of the contour
                    
                    ini_pt = 1 !turn the statement on to satart a new contour for the next segment and element
                    nd3 = first_nd !the third node is the same as the one which open the contour
                    
                ELSE
                
                    x3 = x2 + dx;   y3 = y2 + dy;
                    xnode(n,1) = n; xnode(n,2) = x3; xnode(n,3) = y3;
                    nd3 = n
                    n = n + 1
                    x1 = x3; y1 = y3;
                END IF
                
            ELSE !the segment is not straight
                
                x2 = xc+abs(r)*cos(tet1 + ((e-1)*2.0*dtet + dtet)*sig)
                y2 = yc+abs(r)*sin(tet1 + ((e-1)*2.0*dtet + dtet)*sig)
                xnode(n,1) = n; xnode(n,2) = x2; xnode(n,3) = y2;
                nd2 = n
                n = n + 1
                
                IF (p2 .eq. first_pt .and. e .eq. disc(s,2)) THEN !this elemnt is the last one of the contour
                    
                    ini_pt = 1 !turn the statement on to satart a new contour for the next segment and element
                    nd3 = first_nd !the third node is the same as the one which open the contour
                    
                ELSE
                
                    x3 = xc+abs(r)*cos(tet1 + ((e-1)*2.0*dtet + 2.0*dtet)*sig)
                    y3 = yc+abs(r)*sin(tet1 + ((e-1)*2.0*dtet + 2.0*dtet)*sig)
                    xnode(n,1) = n; xnode(n,2) = x3; xnode(n,3) = y3;
                    nd3 = n
                    n = n + 1
                    x1 = x3;    y1 = y3;
                END IF
                
            END IF
            
            conn(ce,1) = ce;  conn(ce,2) = nd1;  conn(ce,3) = nd2;  conn(ce,4) = nd3; conn(ce,5) = seg(s,5)
            ce = ce + 1
            
        END DO
        
    END DO
    
    !loop over the drops
    do i=1,ndrop
        first_nd = n !first node of the drop
        drop_first_nd = n
        drop_first_last_node(i,1) = i
        drop_first_last_node(i,2) = first_nd
        !loop over the drops elements
        do j=1,drop_nele
            dtet = 2.0*pi/drop_nele
            x1 = drop_geom(i,2) + cos(-dtet*(j-1))*drop_geom(i,4)
            y1 = drop_geom(i,3) + sin(-dtet*(j-1))*drop_geom(i,4)
            x2 = drop_geom(i,2) + cos(-dtet*(j-1)-dtet/2.0)*drop_geom(i,4)
            y2 = drop_geom(i,3) + sin(-dtet*(j-1)-dtet/2.0)*drop_geom(i,4)
            x3 = drop_geom(i,2) + cos(-dtet*(j))*drop_geom(i,4)
            y3 = drop_geom(i,3) + sin(-dtet*(j))*drop_geom(i,4)
            if(j .lt. drop_nele)then
                conn(ce,1) = ce
                conn(ce,2) = n
                xnode(n,1) = n
                xnode(n,2) = x1
                xnode(n,3) = y1
                n = n+1
                conn(ce,3) = n
                xnode(n,1) = n
                xnode(n,2) = x2
                xnode(n,3) = y2
                n = n+1
                conn(ce,4) = n
                conn(ce,5) = 1.0                
            elseif(j .eq. drop_nele)then
                conn(ce,1) = ce
                conn(ce,2) = n
                xnode(n,1) = n
                xnode(n,2) = x1
                xnode(n,3) = y1
                n = n+1
                conn(ce,3) = n
                xnode(n,1) = n
                xnode(n,2) = x2
                xnode(n,3) = y2
                drop_first_last_node(i,3) = n
                n = n+1
                conn(ce,4) = first_nd
                conn(ce,5) = 1.0
            end if
            ce = ce + 1
            
        end do
    end do
    
END SUBROUTINE