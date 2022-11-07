




function Y = cosinorOneFreq(t,theta)
% theta: struct that stores parameters for regression
% t: time of evaluation can be vector
A1=theta.A1; % extract params
phi1=theta.phi1;
f1=theta.f1;
Y=A1.*cos(2*pi*t.*f1-phi1);
end

function [mt_unif,mt_nu] = getSamplingSchedules(NL,NR,tauA,tauB)
Ntimes=NL+NR;
mt1=linspace(0,tauB-tauA,NL+1);
mt1=mt1(1:end-1);
mt2=linspace(tauB-tauA,1,NR+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2]; % construct non-uniform grid 
mt_nu=mod(mt_nu+tauA,1);
mt_unif=linspace(0,1,Ntimes+1); % construct uniform grid 
mt_unif=mt_unif(1:end-1);
end


function Y = cosinorTwoFreq(t,theta)
A1=theta.A1;A2=theta.A2; % extract params
phi1=theta.phi1;phi2=theta.phi2;
f1=theta.f1;f2=theta.f2;
Y=A1.*cos(2*pi*t.*f1-phi1) + A2.*cos(2*pi*t.*T2-phi2); 
end

function X = sampleTruncatedPrior(N,settings)
if nargin<2
	settings.model='cosinorOneFreq';
    settings.run_gpu=false;
	settings.Tprior=100;
    settings.sig1=1;
    settings.sig2=settings.sig1;
    settings.sig3=settings.sig1;
    settings.mu1=1;
    settings.mu2=0;
    settings.mu3=2;
end
T=settings.Tprior;

sig1=settings.sig1;sig2=settings.sig2;sig3=settings.sig3;
mu1=settings.mu1;mu2=settings.mu2;mu3=settings.mu3;

pfun=@(p) (p(:,1)>0).*(p(:,2)>0).*(2*pi>p(:,2)).*(p(:,3)>0).* ...
	normpdf(p,[mu1,mu2,mu3],[sig1,sig2,sig3]); % TODO check this evaluates Gaussian columnwise with correct mean

X=[mu1 mu2 mu3]+zeros(N,3,1); % initial state
Y=randn([N,3,T]).*[sig1 sig2 sig3]; % noise for Brownian motion


if settings.run_gpu
    X=gpuArray(X);
    Y=gpuArray(Y);
end

for t=2:T
    Z=Y(:,:,t)+X;
    diffp=log(pfun(Z))-log(pfun(X));
    logdiff=log(prod(normpdf(X,Z,[sig1 sig2 sig3]),2))- ...
        log(prod(normpdf(Z,X,[sig1 sig2 sig3]),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(N,1));
    X=(rmat<alphamat).*(Z(:,1)>0).*(Z(:,2)>0).*(2*pi>Z(:,2)).*(Z(:,3)>0).*Y(:,:,t) + X;
end

if settings.run_gpu
    X=gather(X);
end

end