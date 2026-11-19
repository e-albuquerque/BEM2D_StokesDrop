function [NOS_2] = rep_no2(No,ELEM,NOS)

    n_el = length(ELEM(:,1));
    NOS_2 = NOS;
    for i=No+1:n_el
       no1 = ELEM(i,2); no2 = ELEM(i,3); no3 = ELEM(i,4);
       x1 = NOS_2(no1,2); y1 = NOS_2(no1,3);
       x2 = NOS_2(no2,2); y2 = NOS_2(no2,3);
       x3 = NOS_2(no3,2); y3 = NOS_2(no3,3);
       
       Xm = (x3+x1)/2; Ym = (y3+y1)/2;
       DX = abs(x3-x1);  DY = abs(y3-y1);
       
       qsi_x = res_eq2(x1,x2,x3,Xm,0);
       qsi_y = res_eq2(y1,y2,y3,Ym,0);
       
       qsi_m = (DX*qsi_x+DY*qsi_y)/(DX+DY);
       
       [N1,N2,N3]=calc_fforma(qsi_m); % Calcula as funções de forma
        N=[N1,0,N2,0,N3,0;
        0,N1,0,N2,0,N3];

        x=N1*x1+N2*x2+N3*x3; % Calcula a coordenada x do ponto de integração
        y=N1*y1+N2*y2+N3*y3; % Calcula a coordenada y do ponto de integração
        
        NOS_2(no2,2) = x ; NOS_2(no2,3) = y ;
        
    end;
  
