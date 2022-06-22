% example settings
clear all 
close all

tps=0:1:24;

settings.osc_fraction=0.05; 
settings.num_samples=1000;
settings.mean_amp=0.4;
settings.noise_to_amp=0.3;
amps=exprnd(settings.mean_amp,1,settings.num_samples);

settings.bimodal_frac=0.9; % 1 means just dist 1, 0 means just dist 2
settings.acro_mu1= pi/2; % phase
settings.acro_mu2= 3*pi/2;
settings.acro_sigma1 = 0.3;
settings.acro_sigma2 = 0.3;
% use CDF to sample bimodal
bimodal_norm_cdf= @(x) settings.bimodal_frac*normcdf(x,settings.acro_mu1,settings.acro_sigma1)...
    +(1-settings.bimodal_frac)*normcdf(x,settings.acro_mu2,settings.acro_sigma2);

rand_samps=rand(1,settings.num_samples);

optim_options = optimset('Display','off');
acro=arrayfun(@(x) fminsearch(@(y) (bimodal_norm_cdf(y) - x)^2,pi,optim_options),rand_samps); %fminbnd works well on bounded domain

subplot(3,1,1)
xlim([0 2*pi])
plot(0:0.01:2*pi,bimodal_norm_cdf(0:0.01:2*pi));

subplot(3,1,2)
hist(acro,sqrt(settings.num_samples))
xlim([0 2*pi])

noise=settings.mean_amp*settings.noise_to_amp*(2*rand(settings.num_samples,length(tps))-1);
%%
tp_mat=ones(settings.num_samples,1)*tps;
acro_mat=acro'*ones(1,length(tps));
phase_mat=tp_mat*2*pi/24-acro_mat;

signal=amps'.*cos(phase_mat)+noise;

subplot(3,1,3)
hold on
for i=1:settings.num_samples
    plot(tps,signal(i,:),'-k','Color',[0 0 0 0.1])
end
hold off


