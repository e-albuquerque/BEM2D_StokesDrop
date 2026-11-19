function [g,b]=calc_gsing(x1,y1,x2,y2,x3,y3,qsi0,mi,qsi,w,...
    gota,lambda,Ca,a_H)

% Pontos de Gauss logarítmico
eta = [0.90425944e-2 0.53971054e-1	0.13531134 0.24705169 0.38021171 ...
    0.52379159 0.66577472 0.79419019 0.89816102 0.96884798];

% Pesos de Gauss Logarítmico
rho = [0.12095474 0.18636310 0.19566066 0.17357723 0.13569597 0.93647084e-1 ...
    0.55787938e-1 0.27159893e-1 0.95151992e-2 0.16381586e-2];

npgr=length(qsi); % Número de pontos de Gauss padrão (regular)
npgs=length(eta); % Número de pontos de Gauss logarítmico (singular)

gs=zeros(2,2); % Inicializa a integral singular
gns=zeros(2,2); % Inicializa a integral não singular
bs=zeros(2,1);
bns=zeros(2,1);

GEL=mi; % Módulo de elasticidade transversal
ni=0.5;
for i=1:npgr % Laço sobre os pontos de integração regular
    [N1,N2,N3]=calc_fforma(qsi(i)); % Funções de forma
    if(qsi0==-1)
        deno=1;
        x0=x1;
        y0=y1;
        ff = [N1 0; % Vetor de funções de forma
            0 N1];
    else
        deno=2;
        x0=x1;
        y0=y1;
        ff = [N2 0; % Vetor de funções de forma
            0 N2];
        
    end
    
    [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(qsi(i)); % Derivadas das funções de forma
    
    dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
    dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);
    
    sx=1/dgamadqsi*dxdqsi; % Component x do vetor tangente
    sy=1/dgamadqsi*dydqsi; % Componente y do vetor tangente
    
    nx=sy; % Componente x do vetor normal
    ny=-sx; % Componente y do vetor normal
    
    [ddN1dqsi,ddN2dqsi,ddN3dqsi]=calc_ddfforma(qsi(i)); % 2as Derivadas das funções de forma
    
    ddxdqsi=ddN1dqsi*x1+ddN2dqsi*x2+ddN3dqsi*x3;
    ddydqsi=ddN1dqsi*y1+ddN2dqsi*y2+ddN3dqsi*y3;
    
    k = abs(dxdqsi*ddydqsi-dydqsi*ddxdqsi)/((dxdqsi^2+dydqsi^2)^(3/2));
    
    x=N1*x1+N2*x2+N3*x3; % Calcula a coordenada x do ponto de integração
    y=N1*y1+N2*y2+N3*y3; % Calcula a coordenada y do ponto de integração
    r=sqrt((x-x0)^2+(y-y0)^2);
    rx=(x-x0)/r;
    ry=(y-y0)/r;
    rd=[rx^2 rx*ry; rx*ry ry^2];
    %dleta = qsi(i)+qsi0
%     l1 = ((qsi(i)+qsi0)*(x1-2*x2+x3)+x3-x1)/deno;
%     l2 = ((qsi(i)+qsi0)*(y1-2*y2+y3)+y3-y1)/deno;
    %modificação Lucas!!!
    l1 = ((qsi(i)+qsi0)*(x1-2*x2+x3)+x3-x1)/deno;
    l2 = ((qsi(i)+qsi0)*(y1-2*y2+y3)+y3-y1)/deno;
    
    lnrns = -log(sqrt(l1^2+l2^2));
    
    % Parte não singular da solução fundamental de deslocamento (u_est_ns[2x2])
    
    u_est_ns = ((3-4*ni)*[lnrns 0 ; 0 lnrns] + rd)/(8*pi*GEL*(1-ni));
    b_est_ns_1 = ([lnrns 0 ; 0 lnrns] + rd)*2/(4*pi*(lambda+1)*Ca);
    %b_est_ns = [b_est_ns_1(1,1)*nx + b_est_ns_1(1,2)*ny
     %               b_est_ns_1(2,1)*nx + b_est_ns_1(2,2)*ny];
    % Iteração de soma da integração não singular
    gns = gns + u_est_ns*ff*dgamadqsi*w(i);
    %bns = bns + b_est_ns*dgamadqsi*w(i);
    bns = bns + b_est_ns_1*[nx ; ny]*k*dgamadqsi*w(i);
end

% Integração singular
for i=1:npgs
    if(qsi0==-1) % O ponto fonte é o primeiro (ou o terceiro, uma vez que
        % a integral é a mesma)
        qsil=2*eta(i)-1;
        [N1,N2,N3]=calc_fforma(qsil); % Funções de forma
        dqsideta=2;
        ff=[N1 0;0 N1];
    else
        qsil=eta(i); % O ponto fonte é o do centro
        [N1,N2,N3]=calc_fforma(qsil); % Funções de forma
        dqsideta=1;
        ff=[N2 0;0 N2];
    end
    [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(qsil); %Derivada das funções de forma
    dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
    dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);
    
    sx=1/dgamadqsi*dxdqsi; % Component x do vetor tangente
    sy=1/dgamadqsi*dydqsi; % Componente y do vetor tangente
    
    nx=sy; % Componente x do vetor normal
    ny=-sx; % Componente y do vetor normal
    
    [ddN1dqsi,ddN2dqsi,ddN3dqsi]=calc_ddfforma(qsil); % 2as Derivadas das funções de forma
    
    ddxdqsi=ddN1dqsi*x1+ddN2dqsi*x2+ddN3dqsi*x3;
    ddydqsi=ddN1dqsi*y1+ddN2dqsi*y2+ddN3dqsi*y3;
    
    k = abs(dxdqsi*ddydqsi-dydqsi*ddxdqsi)/((dxdqsi^2+dydqsi^2)^(3/2));
    
    if(qsi0==-1)
        gs=gs+(3-4*ni)/(8*pi*GEL*(1-ni))*ff*dgamadqsi*rho(i)*dqsideta;
        bs=bs+2/(4*pi*(lambda+1)*Ca)*k*[nx ; ny]*dgamadqsi*rho(i)*dqsideta;
    else
        gs=gs+2*(3-4*ni)/(8*pi*GEL*(1-ni))*ff*dgamadqsi*rho(i)*dqsideta;
        bs=bs+4/(4*pi*(lambda+1)*Ca)*k*[nx ; ny]*dgamadqsi*rho(i)*dqsideta;
    end
end
g=gs+gns;
b=bs+bns;

if (gota==1)
    g = 0.0;
else 
    b = 0.0;
end
