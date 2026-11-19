% Mostrar numeração dos nós
  for no  = 1 : length(PONTOS_INT(:,1))
    n_no = num2str(PONTOS_INT(no,1));
    x_no = PONTOS_INT(no,2);
    y_no = PONTOS_INT(no,3);
    text(x_no,y_no,n_no);
  end;
