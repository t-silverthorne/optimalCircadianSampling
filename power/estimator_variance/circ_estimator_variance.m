clear
close all
addpath('../utils/')
param.useGPU=false;
param.NL=4;
param.NR=4;
param.freq_true=12.3;
param.Amp=1;
param.acro=0;
[t,~]=getSamplingSchedules(param.NL,param.NR,0,0.5);

param.Nreps=1e5;
X=constructX(t,param);
d=param.NL+param.NR;

acro=0:.1:2*pi;
[mu,st]=get_mean_and_std(acro,t,param);
tiledlayout(2,1)
plot(acro,st.amp./mu.amp)
ylim([-1,1])

function [mu,st]=get_mean_and_std(acro,t,param)
d=param.NL+param.NR;
X=constructX(t,param);
Y=param.Amp*cos(2*pi*param.freq_true*t'-acro)+randn(d,numel(acro),param.Nreps,'gpuArray');
beta_hat=pagemtimes((X'*X)\X',Y);

mu.beta=mean(beta_hat,3);
st.beta=std(beta_hat,[],3);

amp=sqrt(abs(sum(beta_hat.^2,1) - beta_hat(1,:,:).^2));
mu.amp=mean(amp,3);
st.amp=std(amp,[],3);% abs val to avoid numerical imag

end