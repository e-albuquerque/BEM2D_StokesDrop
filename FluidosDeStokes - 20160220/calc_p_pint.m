function p=calc_p_pint(ELEM,NOS,PONTOS_INT,desl,trac,mi)
% Calcula a pressão nos pontos internos
k=-1;
n_pint=size(PONTOS_INT,1); % Numero de pontos internos
n_elem=size(ELEM,1); % Numero de elementos
Sx_int=zeros(n_pint,2*n_elem);
Sy_int=zeros(n_pint,2*n_elem);
Dx_int=zeros(n_pint,3*n_elem);
Dy_int=zeros(n_pint,3*n_elem);
tx=zeros(3*n_elem,1);
ty=zeros(3*n_elem,1);
n_pontos_Gauss=8;
[qsi,w]=Gauss_Legendre(-1,1,n_pontos_Gauss);

ind=0;
for i=1:n_elem
    for j=1:3
        ind=ind+1;
        tx(ind)=trac(i,2*j);
        ty(ind)=trac(i,2*j+1);
        el(ind)=i;
        no_el(ind)=j;
    end
end
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
        [dx,dy,sx,sy]=calc_DeS(x1,y1,x2,y2,x3,y3,x_fonte,y_fonte,k,qsi,w); % Chama a
        % functio para cálculo de H e G quando o
        % ponto fonte nao pertence ao elemento
        for nolocal = 1 : 3
            noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
            Sx_int(i,noglobal) = Sx_int(i,noglobal) + sx(nolocal);
            Sy_int(i,noglobal) = Sy_int(i,noglobal) + sy(nolocal);
        end;
        Dx_int(i,3*j-2:3*j) = dx;
        Dy_int(i,3*j-2:3*j) = dy;
    end
end
px=Dx_int*tx-mi*Sx_int*desl(:,2); 
py=Dy_int*ty-mi*Sy_int*desl(:,3); 
p=px+py;