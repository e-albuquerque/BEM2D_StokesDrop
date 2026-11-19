function [NOS_2]=evoluca2(NOS,listu,delta_t,No,ELEM)

    nnos = length(NOS(:,1));
    NOS_2 = NOS;
    for i=(No*2)+1:nnos
       NOS_2(i,2) = NOS_2(i,2)+listu(i*2-1)*delta_t;
       NOS_2(i,3) = NOS_2(i,3)+listu(i*2)*delta_t;
    end 
return