function pout=get_linear_params(pin) 
% a1 = sin coeff period 1
% a2 = cos coeff period 1
% a1 = sin coeff period 2
% a2 = cos coeff period 2
A1=pin(:,1);
A2=pin(:,2);
acro1=pin(:,3);
acro2=pin(:,4);

a1= A1.*cos(acro1);
a2=-A1.*sin(acro1);
a3= A2.*cos(acro2);
a4=-A2.*sin(acro2);

pout=[a1 a2 a3 a4];
end