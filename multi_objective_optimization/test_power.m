clear
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


addpath('utils/')
simtype='medium';
checkUseGPU
param.useGPU=false;
param.useParallel=false;

% waveform
param.NL=4;
param.NR=4;
param.Nmeas=param.NL+param.NR;
param.freq_true=7.4;
param.Amp=2.5;
param.noise=1;

% sample size
param.Nperm=1e2;
param.Nresidual=1e2;
param.Nacro=16; % num. fourier samples
param.method='4tensor';
param.perm_method='fy'; % options fy or randperm

% initial guess
nodes='uniform';
Nmeastot=param.NL+param.NR;
[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

tic
[~,p,~]=simulatePWR_rank4(param,'uniform')
toc


