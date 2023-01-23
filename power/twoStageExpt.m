% periodogram inferred from uniform data, then optimal sampling strategy
% determined based on max of initial periodogram
clear
param.useGPU=false;
param.NL=4;
param.NR=4;
param.freq_true=2; % freq used in regression model
param.Amp=2;
param.noise1=1;
param.noise2=1e-1;

nreps=1e2;
mu=NaN(1,nreps);
for ii=1:nreps
    [~,mu(ii)]=get_min_unif(param);
end
histogram(mu)

function [param,min_unif]=get_min_unif(param)
simtype='rough';
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nresidual=30; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'fast'
        param.Nperm=1e2;
        param.Nresidual=1e2; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'long'
        param.Nperm=1e3;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'verylong'
        param.Nperm=1e3;
        param.Nresidual=5e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
end

acro=0;
[t,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % uniform
Nmeas=param.NL+param.NR;
if param.useGPU
    eps=randn(1,Nmeas,'gpuArray');
else
    eps=randn(1,Nmeas);
end
Y=param.Amp*cos(2*pi*t*param.freq_true-acro)+param.noise1*eps;

[pxx,f]=periodogram(Y,[],[],Nmeas);
[~,mxind]=max(pxx);

param.freq_est=f(mxind);

% estimate uniform power given this period
costfun=@(t) -min(simulatePWR(param,cumsum(abs(t))/sum(abs(t))));
min_unif=min(simulatePWR(param,'uniform'));
end
