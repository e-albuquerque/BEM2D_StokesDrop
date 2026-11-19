function [NOS_2]=evoluca(NOS,listu,delta_t,No,ELEM)

    nnos = length(NOS(:,1));
    NOS_2 = NOS;
    
    n_elem = length(ELEM(:,1));
    
    for i=No+1:n_elem
       no1 = ELEM(i,2); no2 = ELEM(i,3); no3 = ELEM(i,4);
       x1 = NOS(no1,2); y1 = NOS(no1,3);
       x2 = NOS(no2,2); y2 = NOS(no2,3);
       x3 = NOS(no3,2); y3 = NOS(no3,3);
        
      [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(-1); % Calcula as 
                                           % derivadas das funções de forma
    
       dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
       dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
       dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);

       sx=1/dgamadqsi*dxdqsi; % Component x do vetor tangente
       sy=1/dgamadqsi*dydqsi; % Componente y do vetor tangente
    
       nx=sy; % Componente x do vetor normal
       ny=-sx; % Componente y do vetor normal
       
       %nx=-sy; % Componente x do vetor normal
       %ny=sx; % Componente y do vetor normal
        
       prod_esc = [listu(no1*2-1) listu(no1*2)]*[nx ; ny]; 
       
       NOS_2(no1,2) = NOS_2(no1,2)+prod_esc*nx*delta_t;
       NOS_2(no1,3) = NOS_2(no1,3)+prod_esc*ny*delta_t;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(0); % Calcula as 
                                           % derivadas das funções de forma
    
       dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
       dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
       dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);

       sx=1/dgamadqsi*dxdqsi; % Component x do vetor tangente
       sy=1/dgamadqsi*dydqsi; % Componente y do vetor tangente
    
       nx=sy; % Componente x do vetor normal
       ny=-sx; % Componente y do vetor normal
       
       %nx=-sy; % Componente x do vetor normal
       %ny=sx; % Componente y do vetor normal 
       
       prod_esc = [listu(no2*2-1) listu(no2*2)]*[nx ; ny]; 
       
       NOS_2(no2,2) = NOS_2(no2,2)+prod_esc*nx*delta_t;
       NOS_2(no2,3) = NOS_2(no2,3)+prod_esc*ny*delta_t;
    end 
return