%% simulated data
close all
A0=10;
B0=20;
f0=2;
T=500;
Nsamp=1e4;
Nmeas=30;
muA=A0*(1+rand);
sigA=10;
muB=B0*(1+rand);
sigB=10;
sig1=1;sig2=1;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));

% MCMC

X=randn(Nsamp,2);
X(:,1)=X(:,1)*sigA+muA; X(:,2)=X(:,2)*sigB+muB;

pfun=@(par) -sum((y_obs - par(:,1).*cos(2*pi*f0*t_obs) ...
    - par(:,2).*cos(2*pi*f0*t_obs)   ).^2/2,2) - ...
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
close all
tiledlayout(1,2)
nexttile(1)
histogram(X(:,1),80,'normalization','pdf','EdgeColor','none')
xline(A0)
xline(-A0)

nexttile(2)
histogram(X(:,2),80,'normalization','pdf','EdgeColor','none')
xline(B0)
xline(-B0)
%xlim([0 12])

%% variatioanl bayes

