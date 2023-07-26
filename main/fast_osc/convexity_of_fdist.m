addpath('../utils')
Amp=5;
acro=rand(1)*2*pi;
freq=1+rand;
Nmeas=1e2;
alpha=.05;
mt=linspace(0,1,Nmeas);
mu1=zeros(1,Nmeas);
mu1(randsample(1:length(mu1),5,false))=1;
mu2=zeros(1,Nmeas);
mu2(randsample(1:length(mu2),5,false))=1;

F1=getFlambda(Amp,acro,freq,alpha,mt,mu1);
F2=getFlambda(Amp,acro,freq,alpha,mt,mu2);

%%
clf
xv=0:.01:1;
N=8;
x=20;
% plot(xv,ncfpdf(x,2,N-3,xv))
% hold on
% xline(2*x)
% %%
close all
plot(xv,ncfcdf(x,2,N-3,xv),'-k')
%%
hold on
x0=xv(1:end-1);
plot(x0,(-1/2+x0./(2*x0+N-3)).*ncfpdf(x,2,N-3,x0),'--k')
xline(2*x)
yline(0)
%%
N=9
xv=0:.05:3;
x0=xv(1:end-1)
close all
H=hypergeom(0.5*(N-1),1,xv);
plot(xv,H,'-k')
hold on
plot(x0,diff(H),'--k')
%%
xlim([0,10])
%a=(N-1)/2;b=1;
%plot(x0,-0.5*ncfpdf(x,a,b,x0) + )
%%
close all
plot(xv,xv./(2*xv+N-3))
%%
close all
xv=0:.01:Nmeas
Nmeas=9

tiledlayout(2,1)
for ii=0:1:100
    y  = ncfcdf(ii,2,Nmeas-3,xv);
    yp = (y(2:end)-y(1:end-1))/2/.1;
    ypp = (yp(2:end)-yp(1:end-1))/2/.1;
    nexttile(1)
    plot(xv(2:end-1),diff(diff(y))) ;
    hold on
    nexttile(2)
    plot(xv,y) ;
    hold on
    pause(.1)
end

%%
%%
alpha=.05;
Nmeas=8;
lv1=10*rand;
lv2=10*rand;
lv1=4.7;lv2=.5878;gamma=.1636;

LHS=ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,gamma*lv1+(1-gamma)*lv2) ;
RHS= gamma*ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lv1)+...
    (1-gamma)*ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lv2) ;

LHS
RHS

%%
gamma=rand;
gamma*F1+(1-gamma)*F2 > getFlambda(Amp,acro,freq,alpha,mt,gamma*mu1 +(1-gamma)*mu2)
%%
syms t a b
assume(a>b)
int((t^(a-1)/(1-t)^(a-b+1))^2,t,0,1)

function FL = getFlambda(Amp,acro,freq,alpha,mt,mu)
Nmeas  = sum(mu);
csq    = cos(2*pi*freq*mt-acro);
csq    = mu.*csq;
lambda =  Amp^2*csq*csq'; % non-centrality
FL=ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);
end