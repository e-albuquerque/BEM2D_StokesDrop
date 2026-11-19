function [qsi_x] = res_eq2(x1,x2,x3,Xm,qsi_x0)
        qsi_x = 1.0;
    
    while (abs((qsi_x-qsi_x0)/qsi_x) >= 0.001)
        qsi_x0 = qsi_x;
        qsi_x = qsi_x - ((x1-2*x2+x3)*qsi_x^2+(x3-x1)*qsi_x+2*(x2-Xm))...
            /((x1-2*x2+x3)*2*qsi_x+(x3-x1));
    end;

end 