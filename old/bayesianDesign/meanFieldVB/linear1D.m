%% simulated data
clf
rng(12)
A0=10;
B0=7;
f0=1;
T=1000;
Nsamp=1e5;
Nmeas=50;
muA=A0;
sigA=.5;
muB=B0;
sigB=.5;
sig1=.5;sig2=.5;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));

%% MCMC
tic
X=randn(Nsamp,2);
X(:,1)=X(:,1)*sigA+muA; X(:,2)=X(:,2)*sigB+muB;

pfun=@(par) -sum((y_obs - par(:,1).*cos(2*pi*f0*t_obs) ...
    - par(:,2).*sin(2*pi*f0*t_obs)   ).^2/2,2) - ...
    (muA-par(:,1)).^2/2/sigA^2 -  ...
    (muB-par(:,2)).^2/2/sigB^2;

Y=randn(Nsamp,2,T);
Y(:,1)=Y(:,1)*sig1; Y(:,2)=Y(:,2)*sig2;
for t=2:T
    Z=Y(:,:,t)+X;
    diffp=pfun(Z)-pfun(X); 
    logdiff=log(prod(normpdf(X,Z,[sig1 sig2]),2))- ...
        log(prod(normpdf(Z,X,[sig1 sig2]),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(Nsamp,1));
    X=(rmat<alphamat).*Y(:,:,t) + X;
end
toc

%tiledlayout(1,2)
nexttile(1)
histogram(X(:,1),floor(sqrt(Nsamp)),'normalization','pdf','EdgeColor','none')
hold on
xline(A0)

nexttile(2)
histogram(X(:,2),floor(sqrt(Nsamp)),'normalization','pdf','EdgeColor','none')
hold on
xline(B0)

%xlim([0 12])



param.muA=muA;
param.muB=muB;
param.sigA=sigA;
param.sigB=sigB;
param.f0=f0;

% Get maximum of logpdf used for HMC
maxima=fsolve( @(x) wrapLogPDF(x(1),x(2),y_obs,t_obs,param),[1 1]);
for ind=1:2
    nexttile(ind)
    xline(maxima(ind))
end
%% HMC

hmc=hmcSampler(@(x)logpdf(x(1),x(2),y_obs,t_obs,param),[0,0]);
tic
samp=drawSamples(hmc,'NumSamples',Nsamp);
toc


nexttile(1)
histogram(samp(:,1),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none')
nexttile(2)
histogram(samp(:,2),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none')
%% Tuned HMC
hmc=tuneSampler(hmc);
tic
samp=drawSamples(hmc,'NumSamples',Nsamp);
toc


nexttile(1)
histogram(samp(:,1),80,'Normalization','pdf','EdgeColor','none')
nexttile(2)
histogram(samp(:,2),80,'Normalization','pdf','EdgeColor','none')
%% variatioanl bayes
muAINIT=muA;
muBINIT=muB;
sigAINIT=sigA;
sigBINIT=sigB;
tic
muA=1;
muB=1;
sigA=1;
sigB=1;
for i=1:10
    S1A=cos(2*pi*f0*t_obs)*cos(2*pi*f0*t_obs)';
    S2A=sin(2*pi*f0*t_obs)*cos(2*pi*f0*t_obs)';
    S3A=cos(2*pi*f0*t_obs)*y_obs';

    S1B=sin(2*pi*f0*t_obs)*sin(2*pi*f0*t_obs)';
    S2B=sin(2*pi*f0*t_obs)*cos(2*pi*f0*t_obs)';
    S3B=sin(2*pi*f0*t_obs)*y_obs';
    
    muA=(S3A/2 - (S2A*muB)/2 + muAINIT/(2*sigAINIT^2))/(S1A/2 + 1/(2*sigAINIT^2));
    muB=(S3B/2 - (S2B*muA)/2 + muBINIT/(2*sigAINIT^2))/(S1B/2 + 1/(2*sigBINIT^2));
    sigA=1/(S1A + 1/sigAINIT^2)^(1/2);
    sigB=1/(S1B + 1/sigBINIT^2)^(1/2);
end
toc
%%
for ind=1:2
    nexttile(ind)
    xv=0:.1:15;
    if ind==1
        mu=muA;sig=sigA;
    else
        mu=muB;sig=sigB;
    end
    plot(xv,normpdf(xv,mu,sig))
end
% %%
% syms x muA sig S1 S2 S3 muB muINIT sigINIT
% c=fliplr(coeffs(collect((x-muINIT)^2/2/sigINIT^2 + (S1*x^2 +2*x*(muB*S2-S3))/2,x),x));
% m=-c(2)/2/c(1); n=c(3) - (c(2)^2/4/c(1));
% %%
% simplify(c(1)*(x-m)^2+n -( (x-muINIT)^2/2/sigINIT^2 +(S1*x^2 +2*x*(muB*S2-S3))/2  ))
% %%
% sigNew=1/sqrt(2*c(1))
%%
% m

function [p,dp]=logpdf(Avec,Bvec,y_obs,t_obs,param)
muA=param.muA;
muB=param.muB;
sigA=param.sigA;
sigB=param.sigB;
f0=param.f0;
p=-(sum((y_obs - Avec.*cos(2*pi*f0*t_obs) - Bvec.*sin(2*pi*f0*t_obs)   ).^2/2,2) + ...
    (muA-Avec).^2/2/sigA^2 +  ...
    (muB-Bvec).^2/2/sigB^2);
dp(1,:)=-( sum(cos(2.*f0.*t_obs.*pi).*(Avec.*cos(2.*f0.*t_obs.*pi) - y_obs + Bvec.*sin(2.*f0.*t_obs.*pi)),2) + (2.*Avec - 2.*muA)./(2.*sigA.^2));
dp(2,:)=-( sum(sin(2.*f0.*t_obs.*pi).*(Avec.*cos(2.*f0.*t_obs.*pi) - y_obs + Bvec.*sin(2.*f0.*t_obs.*pi)),2) + (2.*Bvec - 2.*muB)./(2.*sigB.^2));
end

function dp = wrapLogPDF(Avec,Bvec,y_obs,t_obs,param)
[~,dp]=logpdf(Avec,Bvec,y_obs,t_obs,param);
end
