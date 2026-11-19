function [g,h,b]=calc_ghnsing(x1,y1,x2,y2,x3,y3,x0,y0,mi,qsi,w...
    ,gota,lambda,Ca,a_H)

% Function que calcula as matrizes g e h locais quando o ponto fonte não
% pertence ao elemento (integração regular)

npg=length(qsi); % Número de pontos de Gauss

h=zeros(2,6); % Inicializa a matriz h do elemento
g=zeros(2,6); % Inicializa a matriz g do elemento
b=zeros(2,1); % Inicializa o vetor b do elemento

for i=1:npg
    [N1,N2,N3]=calc_fforma(qsi(i)); % Calcula as funções de forma
    N=[N1,0,N2,0,N3,0;
        0,N1,0,N2,0,N3];

    x=N1*x1+N2*x2+N3*x3; % Calcula a coordenada x do ponto de integração
    y=N1*y1+N2*y2+N3*y3; % Calcula a coordenada y do ponto de integração
    
    [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(qsi(i)); % Calcula as 
                                           % derivadas das funções de forma
    
    dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
    dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);
    
    [ddN1dqsi,ddN2dqsi,ddN3dqsi]=calc_ddfforma(qsi(i)); % Calcula as 
                                           % derivadas das funções de forma
    ddxdqsi=ddN1dqsi*x1+ddN2dqsi*x2+ddN3dqsi*x3;
    ddydqsi=ddN1dqsi*y1+ddN2dqsi*y2+ddN3dqsi*y3;
    
    sx=1/dgamadqsi*dxdqsi; % Component x do vetor tangente
    sy=1/dgamadqsi*dydqsi; % Componente y do vetor tangente
    
    nx=sy; % Componente x do vetor normal
    ny=-sx; % Componente y do vetor normal
    
    k = abs(dxdqsi*ddydqsi-dydqsi*ddxdqsi)/((dxdqsi^2+dydqsi^2)^(3/2));
    
    [uast,tast,bast]=calc_solfund(x,y,x0,y0,nx,ny,mi,gota...
        ,lambda,Ca,a_H,k); % Calcula as soluções fundamentais
    h=h+tast*N*dgamadqsi*w(i); % Integração da matriz h
    g=g+uast*N*dgamadqsi*w(i); % Integração da matriz g
    b=b+bast*dgamadqsi*w(i); % Integração do vetor b
end
    