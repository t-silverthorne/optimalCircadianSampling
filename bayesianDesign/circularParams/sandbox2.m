%% setup experiment
addpath('utils/')
addpath('FIMs/')
addpath('models/')

NL=3;
NR=8;
tauA=0;
tauB=1/3;
Nparam=100;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

%% get parameter sets
[theta,fnames]=generateParameters(Nparam+1,'cosinorTwoFreq','pseudo-uniform');
ptrue=theta(1,:);
theta=theta(2:end,:);

%% simulate data
Yobs_unif=cosinorTwoFreq(mt_unif,getTheta(ptrue,fnames));
Yobs_nu=cosinorTwoFreq(mt_nu,getTheta(ptrue,fnames));

%% return bayesian parameter estimate

%% find optimal new tauA tauB



function theta=getTheta(tvec,fnames)
tc=num2cell(tvec);
theta=cell2struct(tc(:),fnames);
end









