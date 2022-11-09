% if using Gaussian prior and likelihood, then the posterior is also Gaussian
addpath('../circularParams/utils/')
addpath('../circularParams/FIMs/')
addpath('../circularParams/models/')

%% setup experiment, and get prior from multistart regression

% time grids
NL=5;
NR=5;
tauA=0;
tauB=1/3;
Nparam=100;
settings.dT=0.2;
settings.Ntot=NL+NR;
settings.verbose=false;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);
tobs_mat_nu=[];
Yobs_mat_nu=[];

% get paramsets
model='cosinorOneFreq';
method='pseudo-uniform';
settings.parallel_mode='vectorize';
settings.proposal_method='fixed'; % options: fixed, iterative
settings.FIM_expectation_method='fixed'; % options: fixed, variance
settings.run_gpu=true;
[theta,fnames]=samplePrior(Nparam+1,model,method,settings);
ptrue=theta(1,:);
ptrue(2)=mod(ptrue(2),2*pi); % for more convenient comparison later
if settings.verbose
    ptrue
end

theta=theta(2:end,:);
M=getBayesianFIMcirc(NL+NR,model); % just a function nothing evaluated yet

% simulate measurement
Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));

tobs_mat_nu=mt_nu;
Yobs_mat_nu=Yobs_nu;

% results of multistart regression
bestfit=multiStartRegression(mt_unif,Yobs_unif,model);
[amp_est,acro_est,per_est]=convertToCircularParams(coeffvalues(bestfit),model);

method='ms-prior';
settings.amp_est=amp_est;
settings.acro_est=acro_est;
settings.per_est=per_est;
thetas=samplePrior(100,model,method,settings);
if settings.verbose
    fprintf("%f %f %f\n",amp_est,acro_est,per_est)
    fprintf('%f\n',ptrue)
end
