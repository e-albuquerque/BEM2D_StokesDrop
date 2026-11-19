function [A,b,DESL_PR]  = aplica_CDC(ELEM,NOS,A,B,CDC,matB)
% Troca as colunas das matrizes A e B conforme as condições de contorno

n_el=size(CDC,1);
nnos=size(NOS,1);
todos_valores=zeros(1,6*n_el);
nos_elem=zeros(nnos,4);
cont=0;
DESL_PR=zeros(1,6);

GEL=1;

for el=1:n_el
    no3=ELEM(el,4);
    for no=1:3
        no_global=ELEM(el,no+1);
        if(no==1||no==2)
            nos_elem(no_global,1)=no_global;
            nos_elem(no_global,2)=el;
            nos_elem(no_global,3)=no;
            if(no==1)
                nos_elem(no3,4)=el;
                nos_elem(no3,5)=3;
            end
        end
        for gdl=1:2
            tipoCDC=CDC(el,4*no-3+2*gdl-1);
            valorCDC=CDC(el,4*no-2+2*gdl-1);
            todos_valores(6*el+2*no+gdl-8)=valorCDC;
            if(tipoCDC==0)
                compartilha=0;
                i=1;
                while (~compartilha&&i<=cont&&cont>0)
                    if(no_global==DESL_PR(i,1)&&gdl==DESL_PR(i,2))
                        compartilha=1;
                        DESL_PR(i,5)=el;
                        DESL_PR(i,6)=no;
                    end
                    i=i+1;
                end
                if(~compartilha)
                    cont=cont+1;
                    DESL_PR(cont,1)=no_global;
                    DESL_PR(cont,2)=gdl;
                    DESL_PR(cont,3)=el;
                    DESL_PR(cont,4)=no;
                end
            end
        end
    end
end





n_des_pr = length(DESL_PR(:,1)); % Number of restricted displacements


for de=1 : n_des_pr % for over the restricted nodes
    n_no=DESL_PR(de,1);
    n_gdl=DESL_PR(de,2);
    n_el=DESL_PR(de,3);
    n_no_loc=DESL_PR(de,4);
    ind_A=(2*n_no-2)+n_gdl;
    ind_B=(6*n_el+2*n_no_loc-8)+n_gdl;
    troca = B(:,ind_B);
    B(:,ind_B) = -A(:,ind_A);
    A(:,ind_A) = -troca*GEL;
    if(DESL_PR(de,5)~=0)
        n_el=DESL_PR(de,5);
        n_no_loc=DESL_PR(de,6);
        ind_B=(6*n_el+2*n_no_loc-8)+n_gdl;
        A(:,ind_A)=A(:,ind_A)-GEL*B(:,ind_B);
        B(:,ind_B)=0;
    end;
end;

for i=1:nnos*2;
    vetbB(i,1) = sum(matB(i,:)); 
end 

b=B*todos_valores'+vetbB; % Cálculo do vetor b




