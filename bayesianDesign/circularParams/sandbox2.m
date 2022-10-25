%% setup experiment
clear
addpath('utils/')
addpath('FIMs/')
addpath('models/')

NL=5;
NR=12;
tauA=0;
tauB=1/3;
Nparam=100;
settings.dT=0.2;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

%% get parameter sets
model='cosinorOneFreq';
method='test-spt';
[theta,fnames]=samplePrior(Nparam+1,model,method,settings);
ptrue=theta(1,:);
theta=theta(2:end,:);

M=getBayesianFIMcirc(NL+NR,'cosinorOneFreq');
%%
tic
expectedBayesianFIM(M,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings)
toc
%% simulate data
Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));

% test improper posterior evaluation

evalLikelihood(mt_unif,Yobs_unif,getTheta(ptrue,fnames),model)
evalLikelihood(mt_nu,Yobs_nu,getTheta(ptrue,fnames),model)

%% return bayesian parameter estimate

%% find optimal new tauA tauB


% fminsearch - 
% penalty methods or reparameterization 

% fmincon
% check about rough functions/default optimizer
% BFGS

% global search



%% sandbox
% TODO .1 get rid of
nreps=3;
t_obs_MAT=repmat(mt_nu,nreps,1);
Y_obs_MAT=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(nreps,numel(mt_unif));
%evalLogImproperPosterior(t_obs_MAT,Y_obs_MAT,getTheta(ptrue,fnames),model,method,settings)
thetasamp=samplePosteriorMCMC(1000,fnames,t_obs_MAT,Y_obs_MAT,model,method,settings);

histogram(thetasamp(:,1),80,Normalization="pdf")
hold on
xv=0:.1:30;
plot(xv,exppdf(xv,10))
%plot(xv,exppdf(5,xv))
xline(ptrue(1))
hold off
%histogram(mod(thetasamp(:,2),2*pi))






