syms xi x0 x1 x2

t=transpose(0:0.1:1)

%X=[ones(numel(t),1) cos(2*pi*t) sin(2*pi*t)]
%%
syms theta
assume(theta,'real')
x=@(theta) [1 cos(theta) sin(theta)];

I=@(theta) x(theta)'*x(theta);

I(0)+ I(pi/2)+I(pi) + I(3*pi/2)

%%
n=4;
x= [1 cos(theta) sin(theta)];
xi=arrayfun( @(k) dirac(.1+theta-2*pi*k/n),0:1:(n-1));

int( I(theta)*sum(xi),theta,0,2*pi)