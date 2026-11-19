function mostra_geo(SEGMENTOS,PONTOS,ELEM,NOS)

% Programa de visualização da geometria e malha
nelem=length(ELEM(:,1));
figure
axes('DataAspectRatio',[1,1,1]);
hold on

% Definição da maior dimensão do problema
max_dl = 0;
for lin = 1 : length(SEGMENTOS(:,1))
    p1 = SEGMENTOS(lin,2);
    p2 = SEGMENTOS(lin,3);
    xp1 = PONTOS(p1,2);	yp1 = PONTOS(p1,3);
    xp2 = PONTOS(p2,2);	yp2 = PONTOS(p2,3);
    dl = sqrt((xp1-xp2)^2+(yp1-yp2)^2);
    if dl > max_dl
        max_dl = dl;
    end;
end;


fat = .005*max_dl;		% Fator de escala

% Plota dois pontos em cor branca somente para definir o tamanho da janela
% na qual a geometria será mostrada

emin = min(PONTOS(:,2))-.1*max_dl;
dmax = max(PONTOS(:,2))+.1*max_dl;
imin = min(PONTOS(:,3))-.1*max_dl;
smax = max(PONTOS(:,3))+.1*max_dl;
plot(emin,imin,'w');	% Botton left point
plot(dmax,smax,'w');	% Top right point

% Plota os pontos que definem a geometria como um "x" preto
for i = 1 : length(PONTOS(:,1))
    plot(PONTOS(i,2),PONTOS(i,3),'blackx','markersize',8);
end;

% Gera os pontos sobre os elementos que serão usados na interpolação da
% geometria
ponti = -1 : .2 : 1;
qsi_ponta=[-1,1];
    for p = 1 : 11
        [N1,N2,N3] = calc_fforma(ponti(p));
        FF(p,1)=N1;
        FF(p,2)=N2;
        FF(p,3)=N3;
    end;

for el = 1 : nelem
    
    
    no1 = ELEM(el,2);
    no2 = ELEM(el,3);
    no3 = ELEM(el,4);
    
    
    
    
    x1 = NOS(no1,2);	y1 = NOS(no1,3);
    x2 = NOS(no2,2);	y2 = NOS(no2,3);
    x3 = NOS(no3,2);	y3 = NOS(no3,3);
    
    for i=1:2
        [dN1dqsi,dN2dqsi,dN3dqsi]=calc_dfforma(qsi_ponta(i)); % Calcula as derivadas das
        %    funções de forma
        dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3;
        dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3;
        dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);
        sx=dxdqsi/dgamadqsi; % Component x do vetor tangente na extremidade do elemento
        sy=dydqsi/dgamadqsi; % Componente y do vetor tangente na extremidade do elemento
        n1(i)=sy; % Componente x do vetor normal unitário na extremidade do elemento
        n2(i)=-sx; % Componente y do vetor normal unitário na extremidade do elemento
    end
    
    
    % Compute the PONTOS to interpolate the geometry over the element
    x_el = [x1;x2;x3];
    y_el = [y1;y2;y3];
    XC = FF*x_el;	% Vector (10x1) with the x coordinates of the PONTOS
    YC = FF*y_el;	% Vector (10x1) with the y coordinates of the PONTOS
    
    
    % Plot the element
    plot(XC,YC,'color','black','LineWidth',1.2,'LineStyle','-')
    plot(x_el,y_el,'k.','markersize',6);	% Plot the node of the elements
    % Plot a line in the beginning and in the end of the element
    plot([XC(1)+fat*n1(1) XC(1)-fat*n1(1)],[YC(1)+fat*n2(1) YC(1)-fat*n2(1)],'LineWidth',1.2,'color','black','LineStyle','-');
    plot([XC(end)+fat*n1(2) XC(end)-fat*n1(2)],[YC(end)+fat*n2(2) YC(end)-fat*n2(2)],'LineWidth',1.2,'color','black','LineStyle','-');
end;

