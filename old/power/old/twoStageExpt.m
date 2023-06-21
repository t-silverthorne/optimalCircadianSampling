% periodogram inferred from uniform data, then optimal sampling strategy
% determined based on max of initial periodogram
clear
global simtype
simtype='rough';
checkUseGPU
tic
param.NL=4;
param.NR=4;
param.freq_true=2.5; % freq used in regression model
param.acro_true='rand';
param.Amp=1;
param.noise1=1/3;
param.noise2=1/3;

nreps=5e2;
mu=NaN(1,nreps);
parfor ii=1:nreps
    [~,mu(ii)]=get_min_unif(param);
end
toc
histogram(mu,floor(sqrt(length(mu))))
xlim([0,1])

function [param,min_unif]=get_min_unif(param)
global simtype
checkUseGPU

if isnumeric(param.acro_true)
    acro=param.acro_true;
else
    acro=2*pi*rand;
end

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
