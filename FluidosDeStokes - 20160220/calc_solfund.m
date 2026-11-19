function  [uast,tast,bast]=calc_solfund(x,y,x0,y0,nx,ny,mi,gota...
                                            ,lambda,Ca,a_H,k) 
%Calcula as soluções fundamentais
GEL=mi;
ni=0.5;
R=sqrt((x-x0)^2+(y-y0)^2); % Raio (distância entre ponto fonte e 
                           % ponto campo)
Rx=(x-x0); % Componente x do raio
Ry=(y-y0); % Componente y do raio
% Cálculo de fatores comuns que se repetem nas soluções

rd1=Rx/R;
rd2=Ry/R;
prod1 = 4*pi*(1-ni);
prod12 = 4*pi*(1-ni)*(lambda+1.0)/(1.0-lambda)*(1.0/a_H);
prod2 = (3-4*ni)*log(1/R);
prod3 = rd1*nx + rd2*ny;

if (gota == 0)
    % Solução fundamental de deslocamento
    u11 = (prod2 + rd1^2)/(2*prod1*GEL);
    u22 = (prod2 + rd2^2)/(2*prod1*GEL);
    u12 = (rd1*rd2)/(2*prod1*GEL);

    % Solução fundamental de tração
    p11 = -(prod3*((1-2*ni)+2*rd1^2))/(prod1*R);
    p22 = -(prod3*((1-2*ni)+2*rd2^2))/(prod1*R);
    p12 = -((prod3*2*rd1*rd2)-(1-2*ni)*(rd1*ny-rd2*nx))/(prod1*R);
    p21 = -((prod3*2*rd1*rd2)-(1-2*ni)*(rd2*nx-rd1*ny))/(prod1*R);

    %Solução fundamental para o vetor b da gota
    b1 = 0.0;
    b2 = 0.0;
else
    % Solução fundamental de deslocamento
    u11 = 0.0;
    u22 = 0.0;
    u12 = 0.0;

    % Solução fundamental de tração
    p11 = -(prod3*((1-2*ni)+2*rd1^2))/(prod12*R);
    p22 = -(prod3*((1-2*ni)+2*rd2^2))/(prod12*R);
    p12 = -((prod3*2*rd1*rd2)-(1-2*ni)*(rd1*ny-rd2*nx))/(prod12*R);
    p21 = -((prod3*2*rd1*rd2)-(1-2*ni)*(rd2*nx-rd1*ny))/(prod12*R);

    %Solução fundamental para o vetor b da gota
    b1 = 2*(1/Ca)/(4*pi*(lambda+1))*k*...
       (nx*(prod2 + rd1^2) + ny*(rd1*rd2));
    b2 = 2*(1/Ca)/(4*pi*(lambda+1))*k*...
       (nx*(rd1*rd2) + ny*(prod2 + rd2^2));
end

% Atribuição de uij e pij às variáveis de saída u_est e p_est
uast = [u11 u12
			u12 u22];

tast = [p11 p12
			p21 p22];

bast = [b1
            b2];