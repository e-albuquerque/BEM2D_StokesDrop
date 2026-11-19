clear; % Apaga as variáveis da memória
close all; % Fecha todas as figuras
clc; % Limpa as linhas de comando
dad_2; % Arquivo de entrada de dados

[NOS,ELEM,CDC]=format_dad_parabolico(PONTOS,SEGMENTOS,MALHA,CCSeg); % formata
%  os dados

PONTOS_INT=gera_pontos_internos(NPX,NPY,PONTOS,SEGMENTOS); % gera os
%                   pontos internos

for it=1:10
    
    [G,H,B]=monta_GeH(ELEM,NOS,mi,lambda,Ca,a_H,No); % Monta as matrizes H e G
    
    [matA,b,DESL_PR]=aplica_CDC(ELEM,NOS,H,G,CDC,B); % Aplica as
    %  condições de contorno
    
    x=matA\b; % Calcula as variáveis desconhecidas
    X=x;
    %det = det(matA)
    [listu,listt,u,t] = reordena(x,CDC,DESL_PR); % Separa
    % velocidades de forças de superfície segundo as condições de
    % contorno
    
    u_pint=calc_u_pint(ELEM,NOS,PONTOS_INT,listu,listt,mi,...
        lambda,Ca,a_H,No); %Calcula
    %   as velocidades nos pontos internos
    %p_pint=calc_p_pint(ELEM,NOS,PONTOS_INT,u,t,mi); % Calcula a pressão nos
    % pontos internos
    %map_cor(NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,u_pint) % Faz um mapa de cor
    
    if it==1
        map_cor(NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,u_pint);
    end
    delta_t = 0.5;
    NOS_2 = evoluca2(NOS,listu,delta_t,No,ELEM);
    NOS = NOS_2;
    
    %NOS_2 = evolucao(NOS,listu,delta_t,No,ELEM);
    %NOS = NOS_2;
    
    %NOS_2 = rep_no2(No,ELEM,NOS);
    %NOS = NOS_2;
    
    % NOS_2 = rep_no1(No,ELEM,NOS);
    % NOS = NOS_2;
    
    %NOS_2 = rep_no2(No,ELEM,NOS);
    %NOS = NOS_2;
    map_cor(NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,u_pint);
end
%map_cor(NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,u_pint) % Faz um mapa de cor

% figure
% tri = delaunay(PONTOS_INT(:,2),PONTOS_INT(:,3));
% xi=PONTOS_INT(:,2);
% yi=PONTOS_INT(:,3);
% ci=p_pint;
% for t = 1 : length(tri)	% Percorre todos os triângulos da malha
%
%     x = [xi(tri(t,1));xi(tri(t,2));xi(tri(t,3))];
%     y = [yi(tri(t,1));yi(tri(t,2));yi(tri(t,3))];
%     c = [ci(tri(t,1));ci(tri(t,2));ci(tri(t,3))];
%
%         patch(x,y,c)
%
% end;
%
% shading interp;
