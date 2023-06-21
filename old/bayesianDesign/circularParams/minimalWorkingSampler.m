% plan: fist use MCMC to sample from truncated prior, 
% then use this as initial state of MCMC for sampling from the posterior,
% still ok to use log scaling because you can force the p(x)=0 for x<0 constraint manually 

addpath('utils/')
addpath('FIMs/')
addpath('models/')
%% Data generation step
clear
model='cosinorOneFreq';
fnames={'A1','phi1','T1'};
NL=15; % set up time grid
NR=15;
tauA=0;
tauB=1/3;
settings.dT=0.2;
settings.Ntot=NL+NR;
settings.verbose=false;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

ptrue=[1 0 2];
switch model % simulate measurement
    case 'cosinorOneFreq'
        Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
        %Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));
end       
Yobs_unif_MAT=Yobs_unif;
tobs_unif_MAT=mt_unif;
                tvec=reshape(tobs_unif_MAT,1,numel(tobs_unif_MAT));
                yvec=reshape(Yobs_unif_MAT,1,numel(Yobs_unif_MAT));
%% MCMC step

% GPU
tic
N=1e4; % desired number of samples
T=50; % length of Markov chain
k=3;

mu1=1;mu2=0;mu3=2;
sig1=1;
sig2=sig1;    
sig3=sig1;

logp=@(pmat) -sum((yvec-pmat(:,1).*cos(2*pi*tvec.*pmat(:,3)-pmat(:,2))).^2/2,2) - ...
	(mu1-pmat(:,1)).^2/2/sig1^2 - ...
	(mu2-pmat(:,2)).^2/2/sig2^2 - ...
	(mu3-pmat(:,3)).^2/2/sig3^2 ;

X=sampleTruncatedPrior(N);
Y=randn([N,k,T]);

run_gpu=false;

if run_gpu
    X=gpuArray(X);
    Y=gpuArray(Y);
end

for t=2:T
    Z=Y(:,:,t)+X;
    diffp=logp(Z)-logp(X);
    logdiff=log(prod(normpdf(X,Z,[sig1 sig2 sig3]),2))- ...
        log(prod(normpdf(Z,X,[sig1 sig2 sig3]),2));
    alphamat=min(diffp+logdiff,0);
    rmat=log(rand(N,1));
    X=(rmat<alphamat).*(Z(:,1)>0).*(Z(:,2)>0).*(2*pi>Z(:,2)).*(Z(:,3)>0).*Y(:,:,t) + X;
end
if run_gpu
    X=gather(X);
end
toc
clf

histogram(X(:,1,end),80,'Normalization','pdf')
xv=-8:.01:8;
hold on
%plot(xv,1./(1+(xv/sig).^2)/pi/sig,'linewidth',3)
%ylim([0 0.5])
%xlim([-20 20])
