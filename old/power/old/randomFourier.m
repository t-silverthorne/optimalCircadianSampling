% test if uniform vs non-uniform sampling improves performance when data is
% generated using randomly obtained Fourier coefficients
close all
clear
param.NL=5;
param.NR=3;
param.useGPU=false;
simtype='rough';
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nbatches=30; 
        param.NfourierSamps=50; % num. fourier samples
    case 'fast'
        param.Nperm=1e2;
        param.Nbatches=1e2; 
        param.NfourierSamps=1e2; 
    case 'long'
        param.Nperm=1e3;
        param.Nbatches=3e2; 
        param.NfourierSamps=3e2; 
    case 'verylong'
        param.Nperm=1e3;
        param.Nbatches=1e3; 
        param.NfourierSamps=1e3; 
end

param.per=2; % period used in regression model


param.fourierMeans=[0.9 2]*6;
param.fourierSigma=[.05 .05];
param.NfourierComps=length(param.fourierMeans);

simulatePWR(param,'uniform');
simulatePWR(param,'non-uniform');

function pwr=simulatePWR(param,nodeType)
NL=param.NL;
NR=param.NR;
Nperm=param.Nperm;
Nbatches=param.Nbatches;

NfourierComps=param.NfourierComps;
NfourierSamps=param.NfourierSamps;
fourierMeans=param.fourierMeans;
fourierSigma=param.fourierSigma;

Nmeas=NL+NR;
per=param.per;
switch nodeType
    case 'uniform'
        t=linspace(0,1,Nmeas);
    case 'cheb'
        mc=1:Nmeas;
        t=cos((2*mc-1)*pi/2/Nmeas);
        t=(t+1)/2;
    case 'non-uniform'
        [~,t]=getSamplingSchedules(NL,NR,0,0.25);
        if param.useGPU
            t=gpuArray(t);
        end
end

pwr=0;
for i=1:Nbatches
fourierSigma=reshape(fourierSigma,1,1,length(fourierSigma));
fourierMeans=reshape(fourierMeans,1,1,length(fourierMeans));
if param.useGPU
    permMat=rand(NfourierSamps,Nmeas,Nperm,'gpuArray');
    eps=randn(NfourierSamps,Nmeas,'gpuArray');
    theta=2*pi*rand(NfourierSamps,1,NfourierComps,'gpuArray');
    ampvec=fourierMeans+fourierSigma.*randn(NfourierSamps,1,NfourierComps,'gpuArray');
else
    permMat=rand(NfourierSamps,Nmeas,Nperm);
    eps=randn(NfourierSamps,Nmeas);
    theta=2*pi*rand(NfourierSamps,1,NfourierComps);
    ampvec=fourierMeans+fourierSigma.*randn(NfourierSamps,1,NfourierComps);
end
[~,I]=sort(permMat,2);
pervec=1:NfourierComps;
pervec=reshape(pervec,1,1,length(pervec));

Y=sum(ampvec.*cos(2*pi*t.*pervec-theta),3)+eps;


X=constructX(t,param); % construct linear model

betas_obs=(X'*X)\(X'*Y'); % observed error
fits_obs=(X*betas_obs)';
SSres_obs=sqrt(sum((fits_obs-Y).^2,2));

m=size(Y,1);n=size(Y,2); % generated matrix of permuted samples YI
offMat=repmat((0:m-1)',1,n)*n;
Yp=Y';
YI=pagetranspose(Yp(pagetranspose(I+offMat)));

betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr=pwr+sum(sum(SSres>SSres_obs,3)/Nperm>.95)/NfourierSamps; % pwr is probability of rejecting null hypothesis given theta=theta_0
end
pwr=pwr/Nbatches;
fprintf('%f\n',pwr)
end

function X=constructX(t,param)
per=param.per;
if param.useGPU
    x0=gpuArray(ones(1,length(t)));
    x1=gpuArray(sin(2*pi*per*t));
    x2=gpuArray(cos(2*pi*per*t));
else
    x0=ones(1,length(t));
    x1=sin(2*pi*per*t);
    x2=cos(2*pi*per*t);
end
X= [x0' x1' x2'];

end
