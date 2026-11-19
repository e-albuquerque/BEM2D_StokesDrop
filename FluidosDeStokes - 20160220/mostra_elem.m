% Mostrar numeração dos nós
  for el  = 1 : length(ELEM(:,1))
    n_elem = num2str(ELEM(el,1));
    no2=NOS(ELEM(el,3));
    x_no = NOS(no2,2);
    y_no = NOS(no2,3);
    text(x_no,y_no,n_elem);
  end;
