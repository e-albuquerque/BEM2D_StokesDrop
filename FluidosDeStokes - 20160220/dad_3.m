% Entrada de dados  

% Matriz para definição de pontos 
PONTOS  = [1   0     0
           2   3     0
           3   4     1
           4   6     1
           5   7     0
           6   10    0
           7   10    4
           8   7     4
           9   6     3
           10  4     3
           11  3     4
           12  0     4
           13   3.50  2
           14   4.50  2];
      
% PONTOS  = [1   0   0 
%           2   1   0 
%           3   1   1
%           4   0   1];
          
% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Radio, tipo do elemento] 
% Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
%                          inicial para o ponto final) 
%                   < 0 -> O centro é à direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento é uma linha reta

SEGMENTOS = [1   1   2    0
             2   2   3    0
             3   3   4    0
             4   4   5    0
             5   5   6    0
             6   6   7    0
             7   7   8    0
             8   8   9    0
             9   9   10   0
             10  10  11   0
             11  11  12   0
             12  12  1    0
             13  13  14   -.5
             14  14  13  -.5];

% SEGMENTOS = [1 1 2  0
% 	 	   2 2 3  0
%          3 3 4  0
%          4 4 1  0];

% Matriz para definição da malha

% MALHA = [no do segmento, no de elementos no segmento]

ne=5;
MALHA = [1  ne
	      2  ne 
          3  ne
          4  ne
          5  ne
          6  ne
          7  ne
          8  ne
          9  ne
          10 ne
          11 ne
          12 ne
          13 ne
          14 ne];

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
% Tipo da CDC em x = 3 => (perfil de velocidade parabólico na direção x),
%valor da CDC (velocidade media), distância entre as placas,
%y médio entre as placas.

 CCSeg=[1  0 0 0 0 
        2  0 0 0 0
    	3  0 0 0 0
        4  0 0 0 0
        5  0 0 0 0
        6  1 0 1 0
        7  0 0 0 0
        8  0 0 0 0
        9  0 0 0 0
        10 0 0 0 0
        11 0 0 0 0
        12 3 10 4 2
        13 0 0 0 0
        14 0 0 0 0];
 

mi=1; % Viscosidade do fluido
lambda=1;
Ca = 1.0;
a_H = 2.0*R/Hh;

NPX=11; % Número de pontos internos na direção x
NPY=11; % Número de pontos internos na direção y
