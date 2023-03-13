%
close all
addpath('utils/')
simtype='long';
checkUseGPU
param.useGPU=false;
param.NL=4;
param.NR=4;
param.Nmeas=param.NL+param.NR;
param.freq_true=2.4;
param.Amp=1;
param.noise=.6;
param.Nacro=32;
nodes='uniform';

[acrovec,pwr,est]=wrap_simulatePWR_matperm_fv(param,nodes)

tiledlayout(4,1)
nexttile
plot(acrovec,pwr)
ylim([0,1])
nexttile
plot(acrovec,est.amp_mu+est.amp_st,'-r')
hold on
plot(acrovec,est.amp_mu-est.amp_st,'-r')
ylim([0,4])
yline(param.Amp)
nexttile
plot(acrovec,cos(est.phi_mu))
hold on
plot(acrovec,cos(acrovec))
nexttile
plot(acrovec,est.phi_cvar)
ylim([0,1])
%plot(acrovec,est.amp_mu)
%%
%%% function

%%
addpath('utils')
gpuDevice(1)
simtype='medium';
param.Nperm=1e5;
param.Nresidual=1e1; 
param.useGPU=true;

param.NL=4;
param.NR=4;
param.freq_true=7.6;
param.Amp=1;
param.noise=.5;
param.Nacro=16;

% limitation of my GPU
param.Nperm=1e3;
param.Nresidual=1e3; 
param.Nacro=16;
param.useGPU=true;
[pwr_fast,pwr_medium,pwr_matperm]=run_test_cases(param);

%%
max(abs(pwr_matperm-pwr_slow))

%%
% shows that matrix product is worth it
param.Nperm=1e5;
param.Nresidual=1e1; 
param.Nacro=16;
param.useGPU=true;
run_test_cases(param)

% same samples as before, no GPU (illustrates slowdown)
param.Nperm=1e3;
param.Nresidual=1e3; 
param.Nacro=16;
param.useGPU=false;
run_test_cases(param)

%% test variance
param.Nperm=1e3;
param.Nresidual=1e3; 
param.Nacro=8;
param.useGPU=true;

nn=10;
mvec=NaN(1,nn);
for ii=1:nn
    [~,pwr_matperm]=simulatePWR_matperm_fv(param,'uniform');
    mvec(ii)=min(pwr_matperm);
end
mean(mvec)
std(mvec)

for ii=1:nn
    [~,pwr_matperm]=simulatePWR(param,'uniform');
    mvec(ii)=min(pwr_matperm);
end
mean(mvec)
std(mvec)


function [pwr_medium,pwr_fast,pwr_matperm]=run_test_cases(param)
tic
[~,pwr_medium]=simulatePWR(param,'uniform');
toc
tic
[~,pwr_fast]=simulatePWR_fully_vectorized(param,'uniform');
toc
tic 
[~,pwr_matperm]=simulatePWR_matperm_fv(param,'uniform');
toc
fprintf('\n')
end
%max(abs(pwr_medium-pwr_matperm))
%max(abs(pwr_fast-pwr_matperm))