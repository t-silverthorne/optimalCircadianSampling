% find distance from (x0,y0) to curve f(t) =1/t
close all
clear all
syms x y lambda
assume([x,y],'positive')

x0=2;y0=5;
G=sqrt((x0-x)^2+(y0-y)^2) + lambda*(x*y-1);

S=vpasolve(gradient(G,[x,y,lambda]),[x,y,lambda],[x0 1/x0,1]);
if isempty(S.x)
    S=vpasolve(gradient(G,[x,y,lambda]),[x,y,lambda],[1/y0 y0,1]);
end

% % %%
% % dcurve=@(x0,y0,t) sqrt((x0-t).^2 + (y0-1./t).^2);
% % [num,~]=numden(diff(sqrt((x0-x)^2 + (y0-1/x)^2),x))
% % 
% % F=@(x0,y0) dcurve(x0,y0,roots([1 -x0 0 y0 -1]));
% %%
% close all
% x0=1;
% y0=3;


tvals=0.01:.01:10;
plot(tvals,1./tvals,'-k')
hold on
plot(x0,y0,'ok')
plot(abs(S.x),abs(S.y),'ob')

xlim([0 10])
ylim([0 10])
%%
sym2poly(num,x)