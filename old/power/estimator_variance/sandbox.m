% verify expression for CV is correct
clear
addpath('utils/')
simtype='fast';
checkUseGPU
param.NL=5;
param.NR=3;
param.freq_true=3; % freq used in regression model
param.Amp=1;
param.acro=2;
beta=[0; param.Amp*sin(param.acro); param.Amp*cos(param.acro)];
[t,~]=getSamplingSchedules(param.NL,param.NR,0,0.25);
X=constructX(t,param)

trace(inv(X'*X))
%%
d=length(t);
L=X*((X'*X)\X');
I=eye(d);
sigma=1;

Nreps=1e6;
Tvals=NaN(1,Nreps);

% Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,Nreps);
% beta_hat=(X'*X)\X'*Y;
% CV_MC    = abs(std(beta_hat,[],2)./mean(beta_hat,2));
% 
% CV_exact = abs(sqrt(sigma^2*diag(inv(X'*X)))./beta);
% 
% CV_MC
% CV_exact
acrovec=0:.2:2*pi;
stdvals=NaN(3,numel(acrovec))

for ii=1:numel(acrovec)
    acro=acrovec(ii);
    param.acro=acro;
    beta=[0; param.Amp*sin(param.acro); param.Amp*cos(param.acro)];
    [t,~]=getSamplingSchedules(param.NL,param.NR,0,0.25);
    X=constructX(t,param);
    Y=param.Amp*cos(2*pi*param.freq_true*t'-param.acro)+sigma*randn(d,Nreps);
    beta_hat=(X'*X)\X'*Y;
    stdvals(:,ii)=std(beta_hatX,[],2);
end
close all
plot(acrovec,stdvals)

