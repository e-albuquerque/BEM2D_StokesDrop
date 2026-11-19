function map_cor(NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,desl_pint)


mostra_geo(SEGMENTOS,PONTOS,ELEM,NOS) % mostra a geometria do problema


% hold off
quiver(PONTOS_INT(:,2),PONTOS_INT(:,3),desl_pint(:,2),desl_pint(:,3))




