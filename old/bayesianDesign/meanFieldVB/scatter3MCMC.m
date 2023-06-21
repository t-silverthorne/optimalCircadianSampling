clear
sqrw =@(w) w.^2/sum(w.^2);
sqrw2=@(w) w.^2/sum(w.^2)^2;
A0=2.5;
B0=2.5;
f0=2.1;
fmax=16;
T=1000;
Nsamp=1e4;
Nmeas=8;

t_obs=linspace(0,1,Nmeas+1);
t_obs=t_obs(1:end-1);
y_obs=A0*cos(2*pi*f0*t_obs)+B0*sin(2*pi*f0*t_obs) + randn(1,numel(t_obs));


NgaussA=1;
NgaussB=2;
NgaussT=3;
gradDim=3*(NgaussA+NgaussB+NgaussT);

muAprior  =0;
sigAprior =2;
muBprior  =0;
sigBprior =2;
muTprior  =0;
sigTprior =2;


tic
sig1=2;sig2=2;sig3=2;
X=randn(Nsamp,3);
X(:,1)=X(:,1)*sigAprior+muAprior; X(:,2)=X(:,2)*sigBprior+muBprior;
X(:,1)=X(:,1)*sigTprior+muTprior;

pfun=@(par) -sum((y_obs - par(:,1).*cos(2*pi*0.5*(fmax+fmax*tanh(par(:,3)))*t_obs) ... % likelihood
    - par(:,2).*sin(2*pi*0.5*(fmax +fmax*tanh(par(:,3)))*t_obs)  ).^2/2,2) - ...
    log(sqrt(2*pi)) - ... % normalization of likelihood (sigma=1)
    (muAprior-par(:,1)).^2/2/sigAprior^2 -  ... % prior
    (muBprior-par(:,2)).^2/2/sigBprior^2 - ...  
    (muTprior-par(:,3)).^2/2/sigTprior^2;

Y=randn(Nsamp,3,T);
Y(:,1)=Y(:,1)*sig1; Y(:,2)=Y(:,2)*sig2; Y(:,3)=Y(:,3)*sig3;
for t=2:T
    Z=Y(:,:,t)+X;
    diffp=pfun(Z)-pfun(X); 
    logdiff=log(prod(normpdf(X,Z,[sig1 sig2 sig3]),2))- ...
        log(prod(normpdf(Z,X,[sig1 sig2 sig3]),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(Nsamp,1));
    X=(rmat<alphamat).*Y(:,:,t) + X;
end
toc

scatter3(X(:,1),X(:,2),X(:,3),'.k')