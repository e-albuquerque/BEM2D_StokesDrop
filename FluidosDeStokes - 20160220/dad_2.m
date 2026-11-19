% Entrada de dados  
Hh = 10.0; %altura do canal
Ll = 20.0; %comprimento do canal
R = 2.0; %raio da gota
x_gota = 5.0;
y_gota = 5.0;

% Matriz para definição de pontos 
PONTOS  = [1   0     0
           2   Ll    0
           3   Ll    Hh 
           4   0     Hh
           5   x_gota-R  y_gota
           6   x_gota+R  y_gota];
      
          
% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Radio, tipo do elemento] 
% Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
%                          inicial para o ponto final) 
%                   < 0 -> O centro é à direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento é uma linha reta

SEGMENTOS = [1 1 2  0
             2 2 3  0
             3 3 4  0
             4 4 1  0
             5 5 6  -R
             6 6 5  -R];
n_seg_ext = 4; %numero de segmentos externos
% Matriz para definição da malha

% MALHA = [no do segmento, no de elementos no segmento]

ne=8;
MALHA = [1   2*ne
	      2  ne 
          3  2*ne
          4  ne
          5  ne
          6  ne];

% Condições de contorno nos segmentos
% CCSeg=[Segmento,tipo da CDC em x,valor da CDC em x , ...
%                            tipo da CDC em y,valor da CDC em y]
% tipo da CDC = 0 => a velocidade é conhecida
% tipo da CDC = 1 => a força de superfície é conhecida
% Para condições de contorno de força normal conhecida proceder:
  % tipo da CDC em x = 2, valor da CDC em x = valor da velocidade
  % tangencial
  % Neste caso pode-se atribuir quaisquer valores para tipo da CDC em y e
  % para valor da CDC em y
 CCSeg=[1 0 0 0 0 
        2 1 0 1 0
    	3 0 0 0 0
        4 0 1 0 0
        5 1 0 1 0
        6 1 0 1 0];
 

mi=1; % Viscosidade do fluido
lambda = 1.0;
Ca = 1.0;
a_H = 2.0*R/Hh;
No=0;
Ni=0;

for i=1:length(MALHA(:,1))
    if (i<=n_seg_ext)
        No=No+MALHA(i,2);
    else
        Ni=Ni+MALHA(i,2);
    end  
end


NPX=4; % Número de pontos internos na direção x
NPY=4; % Número de pontos internos na direção y
