function [Dx,Dy,Sx,Sy]=calc_DeS(x1,y1,x2,y2,x3,y3,xd,yd,k,qsi,w)


npg=length(qsi); % Número de pontos de Gauss

Dx=[0 0 0]; % Inicializa a matriz h do elemento
Dy=[0 0 0]; % Inicializa a matriz g do elemento
Sx=[0 0 0]; % Inicializa a matriz h do elemento
Sy=[0 0 0]; % Inicializa a matriz g do elemento


for kk=1:npg
    [N1,N2, N3]=calc_fforma(qsi(kk)); % Calcula as funções de forma
    [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(qsi(kk)); % Calcula as derivadas das
                                                     %    funções de forma    

    x=N1*x1+N2*x2+N3*x3; % Calcula a coordenada x do ponto de integração
    y=N1*y1+N2*y2+N3*y3; % Calcula a coordenada y do ponto de integração
    
    dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
    dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);
    
    sx=dxdqsi/dgamadqsi; % Componente x do vetor tangente
    sy=dydqsi/dgamadqsi; % Componente y do vetor tangente
    nx=sy; % Componente x do vetor normal
    ny=-sx; % Componente y do vetor normal    

    [Tx,Ty,qx,qy]=calc_dsolfund(x,y,xd,yd,nx,ny,k); % Calcula as soluções fundamentais
    Dx=Dx+Tx*[N1,N2, N3]*dgamadqsi*w(kk); % Integral da matriz H
    Dy=Dy+Ty*[N1,N2, N3]*dgamadqsi*w(kk); % Integral da matriz H
    Sx=Sx+qx*[N1,N2, N3]*dgamadqsi*w(kk); % Integral da matriz G
    Sy=Sy+qy*[N1,N2, N3]*dgamadqsi*w(kk); % Integral da matriz G
   
end
return
    