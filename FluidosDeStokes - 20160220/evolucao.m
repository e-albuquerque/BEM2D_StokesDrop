function [NOS_2]=evolucao(NOS,listu,delta_t,No,ELEM)

    nnos = length(NOS(:,1));
    n_elem = length(ELEM(:,1));
    NOS_2 = NOS;
    
    for el = No+1:n_elem
        no1 = ELEM(el,2);
        no2 = ELEM(el,3);
        no3 = ELEM(el,4);
        
        x1 = NOS(no1,2); y1 = NOS(no1,3);
        x2 = NOS(no2,2); y2 = NOS(no2,3);
        x3 = NOS(no3,2); y3 = NOS(no3,3);
        
        %Calculando o vetor normal para o primeiro nó do elemento
        [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(-1);
        
        dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
        dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
        dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);
        
        sx=1/dgamadqsi*dxdqsi; % Component x do vetor tangente
        sy=1/dgamadqsi*dydqsi; % Componente y do vetor tangente
    
        nx=sy; % Componente x do vetor normal
        ny=-sx; % Componente y do vetor normal
        
        udotn = listu(no1*2-1)*nx + listu(no1*2)*ny;
        dxdt = udotn*nx;
        dydt = udotn*ny;
        NOS_2(no1,2) = NOS_2(no1,2)+ dxdt*delta_t;
        NOS_2(no1,3) = NOS_2(no1,3)+ dydt*delta_t;
        
        %Calculando o vetor normal para o segundo nó do elemento
        [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(0);
        
        dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
        dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
        dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);
        
        sx=1/dgamadqsi*dxdqsi; % Component x do vetor tangente
        sy=1/dgamadqsi*dydqsi; % Componente y do vetor tangente
    
        nx=sy; % Componente x do vetor normal
        ny=-sx; % Componente y do vetor normal
        
        udotn = listu(no2*2-1)*nx + listu(no2*2)*ny;
        dxdt = udotn*nx;
        dydt = udotn*ny;
        NOS_2(no2,2) = NOS_2(no2,2)+ dxdt*delta_t;
        NOS_2(no2,3) = NOS_2(no2,3)+ dydt*delta_t;
    end;
   
return