function [G,H,B]=monta_GeH(ELEM,NOS,mi,lambda,Ca,a_H,No)

% Function que monta as matrizes G e H globais

n_el = length(ELEM(:,1));	% Número total de elementos
n_nos = length(NOS(:,1));	% Número total de nós

n_pontos_Gauss=12;
[qsi,w]=Gauss_Legendre(-1,1,n_pontos_Gauss);


% Inicialização das matrizes H e G
H = zeros(2*n_nos,2*n_nos);
G = zeros(2*n_nos,6*n_el);
B = zeros(2*n_nos,n_el);

for i = 1 : n_nos	% Percorre os pontos fontes
    
    % Coordenadas (xf,yf) dos pontos fontes
    x0 = NOS(i,2);
    y0 = NOS(i,3);
    
    for j = 1 : n_el	% Percorre os elementos
        
        % Numeração dos três nós do elemento i
        no1 = ELEM(j,2);
        no2 = ELEM(j,3);
        no3 = ELEM(j,4);
        
        % Coordenadas dos três nós (x1,y1,x2,y2,x3,y3)
        x1 = NOS(no1,2);	y1 = NOS(no1,3);
        x2 = NOS(no2,2);	y2 = NOS(no2,3);
        x3 = NOS(no3,2);	y3 = NOS(no3,3);
        
        % Cálculo das submatrizes h_el e g_el
        
        % Integração não singular
        if (j<=No)
            gota = 0;
            [g,h,b] = calc_ghnsing(x1,y1,x2,y2,x3,y3,x0,y0,mi,qsi,w...
                ,gota,lambda,Ca,a_H);
        else
            gota = 1;
            [g,h,b] = calc_ghnsing(x1,y1,x2,y2,x3,y3,x0,y0,mi,qsi,w...
                ,gota,lambda,Ca,a_H);
        end
        
        % Integração singular
        if ((i == no1) || (i == no3))  % O nó j pertence ao elemento i
            [gsing,bsing]=calc_gsing(x1,y1,x2,y2,x3,y3,-1,mi,qsi,w...
                ,gota,lambda,Ca,a_H); % Integração 
                      %     singular com ponto fonte na extremidade 
                      %     do elemento qsi0=-1 ou qsi0=1)
            if(i==no1)
                g(:,1:2)=gsing;
            else
                g(:,5:6)=gsing;
            end
            b = bsing;
        elseif (i==no2)
            [gsing,bsing]=calc_gsing(x1,y1,x2,y2,x3,y3,0,mi,qsi,w...
                ,gota,lambda,Ca,a_H);% Integração 
                     %     singular com ponto fonte no meio 
                     %     do elemento (qsi0 = 0)
            g(:,3:4)=gsing;
            b = bsing; 
        end
        for nolocal = 1 : 3
            noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
            H(2*i-1:2*i,2*noglobal-1:2*noglobal) =  ...
                H(2*i-1:2*i,2*noglobal-1:2*noglobal) + ...
                h(:,2*nolocal-1:2*nolocal);
        end;
        G(2*i-1:2*i,6*j-5:6*j) = g;
        B(2*i-1:2*i,j) = b;
    end
end

% Calculo dos termos da diagonal da matriz H (consideraçao de corpo a
% temperatura constante)
for m = 1 : n_nos
    H(2*m-1:2*m,2*m-1:2*m) = 0; %zera a diagonal principal
    for n = 1 : n_nos
        if n ~= m
            H(2*m-1:2*m,2*m-1:2*m) = H(2*m-1:2*m,2*m-1:2*m) ...
                - H(2*m-1:2*m,2*n-1:2*n);
        end;
    end;
%     H(2*m-1:2*m,2*m-1:2*m) = H(2*m-1:2*m,2*m-1:2*m) + ...
%     [sum(B(2*m-1,:)) sum(B(2*m-1,:)) ; sum(B(2*m,:)) sum(B(2*m,:))];    
end;

%Lucas!!! Adição dos termos da matriz B nos termos singulares de H
 for m = 1 : n_nos
    H(2*m-1,2*m-1:2*m) = H(2*m-1,2*m-1:2*m) + sum(B(2*m-1,:));
    H(2*m,2*m-1:2*m) = H(2*m,2*m-1:2*m) + sum(B(2*m,:)); 
 end;
return
