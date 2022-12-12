clear
rng(2)


NS=1e3;% number of samples

S=[   -0.5554;
   -0.0348;
   -0.8306];



Ngauss=2;

A0=2.5;
B0=2.5;
f0=2.1;
fmax=15;
Nmeas=12;


x=dlarray(S)
logqlambda(x)

function [y,dydx] = logqlambda(x)
d=3;
muvec=dlarray([0.2272    0.2069;
    0.6967    0.3203;
    0.3260    0.4498]);

Lvec =dlarray([0.0346    0.7466;
    0.6136   -0.0204;
    0.2844    0.4052;
    0.6853   -0.0347;
    0.9540    0.5281;
    0.3942   -0.0035]);

wvec=dlarray([0.5987
    0.4013]);

Sigm1=ivh(Lvec(:,1))*ivh(Lvec(:,1))';
Sigm2=ivh(Lvec(:,2))*ivh(Lvec(:,2))';

dS1=det(Sigm1);
dS2=det(Sigm2);
invSigm1=dlarray(inv(Sigm1));
invSigm2=dlarray(inv(Sigm2));
Sigm1=dlarray(Sigm1);
Sigm2=dlarray(Sigm2);



y=log(wvec(1).^2/sum(wvec.^2)*(sqrt(dS1)^(-1)*(2*pi)^(-d/2)*...
               exp(-0.5*(x-muvec(:,1))'*(invSigm1*(x-muvec(:,1)) ))) + ...
wvec(2).^2/sum(wvec.^2)*(sqrt(dS2)^(-1)*(2*pi)^(-d/2)*...
               exp(-0.5*(x-muvec(:,2))'*(invSigm2*(x-muvec(:,2))))))
dydx=dlgradient(y,x);
end

function Lmat = ivh(Lvec)
% inverse of vech map
d= (8*length(Lvec) + 1)^(1/2)/2 - 1/2;
Lmat = zeros(d,d);
i0=1;
for i=1:d
    Lmat(i:d,i)=Lvec(i0:i0+d-i);
    i0=i0+d+1-i;
end
end