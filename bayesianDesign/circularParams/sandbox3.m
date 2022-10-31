clear
addpath('utils/')
addpath('FIMs/')
addpath('models/')

%% setup experiment

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
%% first experiment, generate sample from prior
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

%%
method='test-spt';
settings.speed='slow';% options: slow, fast
settings.run_gpu=false;
tic
samplePosteriorMCMC(1e4,fnames,tobs_mat_nu, ...
                Yobs_mat_nu,model,method,settings);
toc 
%% optimization step
costfunFIM = @(tau,delta) wrapExpectedBayesianFIM(NL,NR,tau,tau+delta,M,fnames,tobs_mat_nu, ...
                                Yobs_mat_nu,model,method,settings);

%%

opts = optimoptions(@particleswarm,'HybridFcn',@fminsearch,...
                                   'Display','iter', ...
                                   'SwarmSize',300, ...
                                   'UseParallel',true); %TODO: set true
sampstrat=particleswarm(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,fnames,tobs_mat_nu, ...
                                Yobs_mat_nu,model,method,settings), ...
                                2,[],[],opts);


%costfunFIM(0,0.1)
%costfunFIM(0,1/2)


function C=wrapExpectedBayesianFIM(NL,NR,tauA,tauB,M,fnames,tobs_mat, ...
                                Yobs_mat,model,method,settings)
tauA=1+0.5*tanh(tauA);
tauB=1+0.5*tanh(tauB);
[~,tmeas_prop]=getSamplingSchedules(NL,NR,tauA,tauB);
C=expectedBayesianFIM(M,fnames,tmeas_prop,tobs_mat,Yobs_mat,model,method,settings);
end


%%

function thetanew=remapTheta(theta,model)
switch model
    case 'cosinorOneFreq'
        thetanew(1)=theta(1).^2;
        thetanew(2)=pi*(1+tanh(theta(2)));
        thetanew(3)=0.5*(1+tanh(theta(3)));
end
end


% ptrue
% remapTheta(thetaest,model)
% settings.cint=confint(bestfit);
% settings.coeffs=coeffvalues(bestfit);
% diff(settings.cint,1)/2
% 
% 
% function error=costFunPS(theta,t_obs,Y_obs,fnames,model)
% theta(1)=theta(1).^2;
% theta(2)=pi*(1+tanh(theta(2)));
% theta(3)=0.5*(1+tanh(theta(3)));
% error= sum((Y_obs - model(t_obs,getTheta(theta,fnames))).^2);
% end
