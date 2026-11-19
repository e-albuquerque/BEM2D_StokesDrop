!This subrouitne relocate the node just like it is exposed in Wrobel 2009
    !Drop deformation in Stokes flow through converging channels
    
subroutine nodes_relocation
    use global
    implicit none
    integer :: el, nd, drp
    integer :: num_drp, ii
    real :: x1(2), x2(2), x3(2)
    real :: x1r(2), x2r(2), x3r(2), x1l(2), x2l(2), x3l(2)
    real :: Xm, Ym, dX, dY
    real :: a1, b1, c1, a2, b2, c2
    real :: qsix, qsiy, qsix_old, qsiy_old, qsi_m
    real :: err
    real :: x2n(2)
    real :: Tol
    
    Tol = 10.0**(-3)
    
    !relocating the second local nodes
    do el = drop_f_elem,nele
        x1(:) = xnode(conn(el,2),2:3)
        x2(:) = xnode(conn(el,3),2:3)
        x3(:) = xnode(conn(el,4),2:3)
        
        Xm = (x3(1)+x1(1))/2.0
        Ym = (x3(2)+x1(2))/2.0
        dX = abs(x3(1)-x1(1))
        dY = abs(x3(2)-x1(2))
        
        a1 = x1(1)-2.0*x2(1)+x3(1)
        b1 = x3(1)-x1(1)
        c1 = 2.0*(x2(1)-Xm)
        
        a2 = x1(2)-2.0*x2(2)+x3(2)
        b2 = x3(2)-x1(2)
        c2 = 2.0*(x2(2)-Ym)
        
        qsix = 0.0d0
        qsiy = 0.0d0
        err = 100.0
        do while(err .gt. Tol)
            qsix_old = qsix;
            
            qsix = qsix_old - (a1*qsix**2 + b1*qsix + c1)/(2.0*a1*qsix + b1)
            
            err = abs(qsix-qsix_old)    
        end do
        
        err = 100.0
        do while(err .gt. Tol)
            qsiy_old = qsiy;
            
            qsiy = qsiy_old - (a2*qsiy**2 + b2*qsiy + c2)/(2.0*a2*qsiy + b2)
            
            err = abs(qsiy-qsiy_old)
        end do
        
        qsi_m = (dX*qsix + dY*qsiy)/(dX + dY)
        
        x2n(1) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(1) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(1) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(1)
        x2n(2) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(2) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(2) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(2)
        
        if(qsi_m .lt. 0.9 .or. qsi_m .gt. -0.9)then
            xnode(conn(el,3),2:3) = x2n
        end if
    end do
    
    !relocating the third local nodes
    do el = drop_f_elem,nele
        x1(:) = xnode(conn(el,3),2:3)
        x2(:) = xnode(conn(el,4),2:3)
        
        num_drp = 0
        do drp=1,ndrop
            if(conn(el,4) .eq. drop_first_last_node(drp,2))then
                num_drp = drp
            end if
        end do
        
        if(num_drp .eq. 0)then
            x3(1:2) = xnode(conn(el+1,3),2:3)
        else 
            x3(1:2) = xnode(drop_first_last_node(num_drp,2)+1,2:3)
        end if
        
        Xm = (x3(1)+x1(1))/2.0
        Ym = (x3(2)+x1(2))/2.0
        dX = abs(x3(1)-x1(1))
        dY = abs(x3(2)-x1(2))
        
        a1 = x1(1)-2.0*x2(1)+x3(1)
        b1 = x3(1)-x1(1)
        c1 = 2.0*(x2(1)-Xm)
        
        a2 = x1(2)-2.0*x2(2)+x3(2)
        b2 = x3(2)-x1(2)
        c2 = 2.0*(x2(2)-Ym)
        
        qsix = 0.0d0
        qsiy = 0.0d0
        err = 100.0
        do while(err .gt. Tol)
            qsix_old = qsix;
            
            qsix = qsix_old - (a1*qsix**2 + b1*qsix + c1)/(2.0*a1*qsix + b1)
            
            err = abs(qsix-qsix_old)
        end do
        
        err = 100.0
        do while(err .gt. Tol)
            qsiy_old = qsiy;
            
            qsiy = qsiy_old - (a2*qsiy**2 + b2*qsiy + c2)/(2.0*a2*qsiy + b2)
            
            err = abs(qsiy-qsiy_old)
        end do
        
        qsi_m = (dX*qsix + dY*qsiy)/(dX + dY)
        
        x2n(1) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(1) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(1) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(1)
        x2n(2) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(2) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(2) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(2)
        
        if(qsi_m .lt. 0.9 .or. qsi_m .gt. -0.9)then
            xnode(conn(el,4),2:3) = x2n
        end if   
    end do
    
    ! do el = drop_f_elem,nele
    !    x1r(:) = xnode(conn(el,3),2:3)
    !    x2r(:) = xnode(conn(el,4),2:3)
    !    
    !    num_drp = 0
    !    do drp=1,ndrop
    !        if(conn(el,4) .eq. drop_first_last_node(drp,2))then
    !            num_drp = drp
    !        end if
    !    end do
    !    
    !    if(num_drp .eq. 0)then
    !        x3r(1:2) = xnode(conn(el+1,3),2:3)
    !    else 
    !        x3r(1:2) = xnode(drop_first_last_node(num_drp,2)+1,2:3)
    !    end if
    !    
    !    x3l(:) = xnode(conn(el,3),2:3)
    !    x2l(:) = xnode(conn(el,2),2:3)
    !    
    !    num_drp = 0
    !    do drp=1,ndrop
    !        if(conn(el,2) .eq. drop_first_last_node(drp,2))then
    !            num_drp = drp
    !        end if
    !    end do
    !    
    !    if(num_drp .eq. 0)then
    !        x1l(1:2) = xnode(conn(el-1,3),2:3)
    !    else 
    !        x1l(1:2) = xnode(drop_first_last_node(num_drp,3),2:3)
    !    end if
    !    
    !    Xm = (x3(1)+x1(1))/2.0
    !    Ym = (x3(2)+x1(2))/2.0
    !    dX = abs(x3(1)-x1(1))
    !    dY = abs(x3(2)-x1(2))
    !    
    !    a1 = x1(1)-2.0*x2(1)+x3(1)
    !    b1 = x3(1)-x1(1)
    !    c1 = 2.0*(x2(1)-Xm)
    !    
    !    a2 = x1(2)-2.0*x2(2)+x3(2)
    !    b2 = x3(2)-x1(2)
    !    c2 = 2.0*(x2(2)-Ym)
    !    
    !    qsix = 0.0d0
    !    qsiy = 0.0d0
    !    err = 100.0
    !    do while(err .gt. Tol)
    !        qsix_old = qsix;
    !        
    !        qsix = qsix_old - (a1*qsix**2 + b1*qsix + c1)/(2.0*a1*qsix + b1)
    !        
    !        err = abs(qsix-qsix_old)
    !    end do
    !    
    !    err = 100.0
    !    do while(err .gt. Tol)
    !        qsiy_old = qsiy;
    !        
    !        qsiy = qsiy_old - (a2*qsiy**2 + b2*qsiy + c2)/(2.0*a2*qsiy + b2)
    !        
    !        err = abs(qsiy-qsiy_old)
    !    end do
    !    
    !    qsi_m = (dX*qsix + dY*qsiy)/(dX + dY)
    !    
    !    x2n(1) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(1) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(1) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(1)
    !    x2n(2) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(2) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(2) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(2)
    !    
    !    xnode(conn(el,4),2:3) = x2n
    !end do
    
    !relocating the second local nodes
    do el = drop_f_elem,nele
        x1(:) = xnode(conn(el,2),2:3)
        x2(:) = xnode(conn(el,3),2:3)
        x3(:) = xnode(conn(el,4),2:3)
        
        Xm = (x3(1)+x1(1))/2.0
        Ym = (x3(2)+x1(2))/2.0
        dX = abs(x3(1)-x1(1))
        dY = abs(x3(2)-x1(2))
        
        a1 = x1(1)-2.0*x2(1)+x3(1)
        b1 = x3(1)-x1(1)
        c1 = 2.0*(x2(1)-Xm)
        
        a2 = x1(2)-2.0*x2(2)+x3(2)
        b2 = x3(2)-x1(2)
        c2 = 2.0*(x2(2)-Ym)
        
        qsix = 0.0d0
        qsiy = 0.0d0
        err = 100.0
        do while(err .gt. Tol)
            qsix_old = qsix;
            
            qsix = qsix_old - (a1*qsix**2 + b1*qsix + c1)/(2.0*a1*qsix + b1)
            
            err = abs(qsix-qsix_old)
        end do
        
        err = 100.0
        do while(err .gt. Tol)
            qsiy_old = qsiy;
            
            qsiy = qsiy_old - (a2*qsiy**2 + b2*qsiy + c2)/(2.0*a2*qsiy + b2)
            
            err = abs(qsiy-qsiy_old)
        end do
        
        qsi_m = (dX*qsix + dY*qsiy)/(dX + dY)
        
        x2n(1) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(1) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(1) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(1)
        x2n(2) = 0.5d0*qsi_m*(qsi_m-1.0d0)*x1(2) + (1.0d0-qsi_m)*(1.0d0+qsi_m)*x2(2) + 0.5d0*qsi_m*(qsi_m+1.0d0)*x3(2)
        
        if(qsi_m .lt. 0.9 .or. qsi_m .gt. -0.9)then
            xnode(conn(el,3),2:3) = x2n
        end if   
        
    end do
    
end subroutine nodes_relocation