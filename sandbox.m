% example settings
clear all 
close all

tps=0:1:24;
% stngs short for settings
stngs.osc_fraction=0.05; 
stngs.num_samples=100;
stngs.mean_amp=0.4;
stngs.noise_to_amp=0.3;
amps=exprnd(stngs.mean_amp,1,stngs.num_samples);

stngs.bimodal_frac=0.9; % 1 means just dist 1 0 means just dist 2
stngs.acro_mu1= pi/2; % phase
stngs.acro_mu2= 3*pi/2;
stngs.acro_sigma1 = 0.3;
stngs.acro_sigma2 = 0.3;
% use CDF to sample bimodal
bimodal_norm_cdf= @(x) stngs.bimodal_frac*normcdf(x,stngs.acro_mu1,stngs.acro_sigma1)...
    +(1-stngs.bimodal_frac)*normcdf(x,stngs.acro_mu2,stngs.acro_sigma2);

rand_samps=rand(1,stngs.num_samples);

optim_options = optimset('Display','off');
acro=arrayfun(@(x) fminsearch(@(y) abs(bimodal_norm_cdf(y) - x),pi,optim_options),rand_samps);

subplot(3,1,1)
xlim([0 2*pi])
plot(0:0.01:2*pi,bimodal_norm_cdf(0:0.01:2*pi));
subplot(3,1,2)
hist(acro,sqrt(stngs.num_samples))
xlim([0 2*pi])

noise=stngs.mean_amp*stngs.noise_to_amp*(rand(stngs.num_samples,length(tps))-0.5);

tp_mat=ones(stngs.num_samples,1)*tps;
acro_mat=acro'*ones(1,length(tps));
phase_mat=tp_mat*2*pi/24-acro_mat;

signal=amps'.*cos(phase_mat)+noise;

subplot(3,1,3)
hold on
for i=1:stngs.num_samples
    plot(tps,signal(i,:),'-k','Color',[0 0 0 0.1])
end
hold off





%%
plot(tps,stngs.amp*cos(tps*2*pi/24)+stngs.amp*stngs.noise_to_amp*(rand(1,length(tps))-0.5),'ok')
hold on
plot(0:0.1:24,stngs.amp*cos((0:0.1:24)*2*pi/24),'-b')
hold off
%%