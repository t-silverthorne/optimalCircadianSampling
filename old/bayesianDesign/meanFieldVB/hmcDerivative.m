A0=5;
B0=2;
f0=1/2;
T=5000;
Nsamp=1e4;
Nmeas=10;
muA=A0*(1+.1*rand);
sigA=1;
muB=B0*(1+.1*rand);
sigB=1;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));

%% HMC
param.muA=muA;
param.muB=muB;
param.sigA=sigA;
param.sigB=sigB;
param.f0=f0;
logpdf([0;1],[1;0],y_obs,t_obs,param)

hmc=hmcSampler(@(x)logpdf(x(1),x(2),y_obs,t_obs,param),[0,0]);
tic
samp=drawSamples(hmc,'NumSamples',Nsamp);
toc
%%


%%

nexttile(1)
histogram(samp(:,1),80,'Normalization','pdf','EdgeColor','none')
nexttile(2)
histogram(samp(:,2),80,'Normalization','pdf','EdgeColor','none')
%%

syms f0 t_obs y_obs Avec Bvec muA muB sigA sigB

%(y_obs - Avec.*cos(2*pi*f0*t_obs) - Bvec.*sin(2*pi*f0*t_obs) + (muA-Avec).^2/2/sigA^2 + (muB-Bvec).^2/2/sigB^2

dpdA1=vectorize(diff( (y_obs - Avec.*cos(2*pi*f0*t_obs) - Bvec.*sin(2*pi*f0*t_obs)   ).^2/2  ,Avec))
dpdA2=vectorize(diff( (muA-Avec).^2/2/sigA^2 + (muB-Bvec).^2/2/sigB^2 ,Avec))

dpdB1=vectorize(diff( (y_obs - Avec.*cos(2*pi*f0*t_obs) - Bvec.*sin(2*pi*f0*t_obs)   ).^2/2  ,Bvec))
dpdB2=vectorize(diff( (muA-Avec).^2/2/sigA^2 + (muB-Bvec).^2/2/sigB^2 ,Bvec))


 


function [p,dp]=logpdf(Avec,Bvec,y_obs,t_obs,param)
muA=param.muA;
muB=param.muB;
sigA=param.sigA;
sigB=param.sigB;
f0=param.f0;
p=-sum((y_obs - Avec.*cos(2*pi*f0*t_obs) - ...
    Bvec.*cos(2*pi*f0*t_obs)   ).^2/2,2) - ...
    (muA-Avec).^2/2/sigA^2 -  ...
    (muB-Bvec).^2/2/sigB^2;
dp(1,:)=sum(-cos(2.*f0.*t_obs.*pi).*(Avec.*cos(2.*f0.*t_obs.*pi) - ...
    y_obs + Bvec.*cos(2.*f0.*t_obs.*pi)),2) -(2.*Avec - 2.*muA)./(2.*sigA.^2); 
dp(2,:)=sum(-cos(2.*f0.*t_obs.*pi).*(Avec.*cos(2.*f0.*t_obs.*pi) - ...
    y_obs + Bvec.*cos(2.*f0.*t_obs.*pi)),2) -(2.*Bvec - 2.*muB)./(2.*sigB.^2);
end
