x=dlarray(-10);
[~,dydx1]=dlfeval(@wgauss,x);
[~,dydx2]=dlfeval(@wgauss2,x);
[~,dydx3]=dlfeval(@wgauss3,x);
fprintf('\n')
fprintf('%d\n',extractdata(dydx1))
fprintf('%d\n',extractdata(dydx2))
fprintf('%d\n',extractdata(dydx3))

function [y,dydx] = wgauss(x)
c1=0.2;
c2=0.8;
a1=10;
a2=5;
y= log(c1*exp(-a1*x^2) + c2*exp(-a2*x^2));
dydx=dlgradient(y,x);
end

function [y,dydx] = wgauss2(x)
c1=0.2;
c2=0.8;
a1=10;
a2=5;
y= log(c1*exp(a2*x^2) + c2*exp(a1*x^2))-a1*x^2 -a2*x^2;
dydx=dlgradient(y,x);
end


function [y,dydx] = wgauss3(x)
c1=0.2;
c2=0.8;
a1=10;
a2=5;
y= log(1/(c1*exp(-a1*x^2) + c2*exp(-a2*x^2)));
dydx=-dlgradient(y,x);
end