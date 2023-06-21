%syms x1 x2
close all
clear
f=@(x1,x2) sqrt(1./x1.^2 +1./x2.^2)./(1./x1+1./x2)

x=0:.05:4;
[X,Y]=meshgrid(x,x)

surf(X,Y,f(X,Y))
%%
syms x1 x2
f=sqrt(1/x1.^2 +1/x2.^2)./(1/x1+1/x2)
simplify(gradient(f))

