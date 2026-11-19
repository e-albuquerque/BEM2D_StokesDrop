function desl_int=calc_desl_int(ELEM,NOS,PONTOS_INT,desl,trac,mi)

% Calcula os deslocamentos nos pontos internos


n_pontos_Gauss=8;
[qsi,w]=Gauss_Legendre(-1,1,n_pontos_Gauss);


n_pint=length(PONTOS_INT(:,1)); % Numero de pontos internos
n_elem=length(ELEM(:,1)); % Numero de elementos
H_int=zeros(2*n_pint,4*n_elem);
G_int=zeros(2*n_pint,6*n_elem);
desl_int=zeros(n_pint,1);
for i=1:n_pint % Laço sobre os pontos internos
    x_fonte=PONTOS_INT(i,2); % Coordenada x do ponto fonte
    y_fonte=PONTOS_INT(i,3); % Coordenada y do ponto fonte
    for j=1:n_elem
        x1=NOS(ELEM(j,2),2); % Coordenada x do inicio do elemento
        y1=NOS(ELEM(j,2),3); % Coordenada y do inicio do elemento
        x2=NOS(ELEM(j,3),2); % Coordenada x do final do elemento
        y2=NOS(ELEM(j,3),3); % Coordenada y do final do elemento
        x3=NOS(ELEM(j,4),2); % Coordenada x do final do elemento
        y3=NOS(ELEM(j,4),3); % Coordenada y do final do elemento
        [g,h] = calc_ghnsing(x1,y1,x2,y2,x3,y3,x_fonte,y_fonte,mi, ...
            qsi,w);%Chama a function para cálculo de H e G quando o
        % ponto fonte nao pertence ao elemento
        for nolocal = 1 : 3
            noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
             H_int(2*i-1:2*i,2*noglobal-1:2*noglobal) =  ...
                H_int(2*i-1:2*i,2*noglobal-1:2*noglobal) + ...
                h(:,2*nolocal-1:2*nolocal);
        end;
        G_int(2*i-1:2*i,6*j-5:6*j) = g;
    end
end

vetdesl_pint=G_int*trac-H_int*desl; % Vetor que contem os deslocamentos 
                      % nos pontos internos
for i=1:n_pint % Laço sobre os pontos internos
   desl_int(i,1)=1;
   desl_int(i,2)=vetdesl_pint(2*i-1);
   desl_int(i,3)=vetdesl_pint(2*i);
end
return
   
                       