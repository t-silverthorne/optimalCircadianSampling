clear
addpath('utils/') 
addpath('FIMs/') 
addpath('models/')

%% get estimate of parameters from uniform experiment
NL=5; % set up time grid
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
M=getBayesianFIMcirc(NL+NR,model); % just a function nothing evaluated yet

model='cosinorOneFreq'; % get true parameter value
method='pseudo-uniform';
[theta,fnames]=samplePrior(Nparam+1,model,method,settings);
ptrue=theta(1,:);
ptrue(2)=mod(ptrue(2),2*pi); 
if settings.verbose
    ptrue
end
theta=theta(2:end,:);

switch model % simulate measurement
    case 'cosinorOneFreq'
        Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
        %Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));
end       

% use results of multistart regression to construct prior
bestfit=multiStartRegression(mt_unif,Yobs_unif,model);
[amp_est,acro_est,per_est]=convertToCircularParams(coeffvalues(bestfit),model);
settings.amp_est=amp_est;
settings.acro_est=acro_est;
settings.per_est=per_est;

settings.speed='fast';% options: slow, fast
settings.run_gpu=false;
settings.parallel_mode='vectorize';
settings.prop='fixed';
settings.prop   
settings.FIM_expectation_method='variance';
settings.batch_size=1e4;
settings.var_cut=1e-1;

%% particle swarm
opts = optimoptions(@particleswarm,'HybridFcn',@fminsearch,...
                                   'Display','iter', ...
                                   'SwarmSize',100, ...
                                   'UseParallel',true, ...
                                   'PlotFcn','pswplotbestf'); %TODO: set true
%%

sampstrat=particleswarm(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,fnames,tobs_mat_nu, ...
                                Yobs_mat_nu,model,method,settings), ...
                                2,[],[],opts);
%% fminsearch
opts=optimset('PlotFcns',@optimplotfval,'MaxFunEvals',100,'UseParallel',true)
fminsearch(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,fnames,tobs_mat_nu, ...
                                Yobs_mat_nu,model,method,settings),[0 1/3],opts)
%% simulanneal
opts=optimoptions(@simulannealbnd,'Display','iter', ...
                    'FunctionTolerance',1e-2);
sampstrat=simulannealbnd(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,fnames,tobs_mat_nu, ...
                                Yobs_mat_nu,model,method,settings),[0 1/2]);
%% fmincon
opts=optimset('PlotFcns',@optimplotfval)
fmincon(@(x) -wrapExpectedBayesianFIMConstrained(NL,NR,x(1),x(1)+x(2),M,fnames,tobs_mat_nu, ...
                                Yobs_mat_nu,model,method,settings), ...
                                [0 1/2],[],[],[],[],[0 0],[1 1], ...
                                 [],opts)

function C=wrapExpectedBayesianFIM(NL,NR,tauA,tauB,M,fnames,tobs_mat, ...
                                Yobs_mat,model,method,settings)
tauA=1+0.5*tanh(tauA);
tauB=1+0.5*tanh(tauB);
[~,tmeas_prop]=getSamplingSchedules(NL,NR,tauA,tauB);
C=expectedBayesianFIM(M,fnames,tmeas_prop,tobs_mat,Yobs_mat,model,method,settings);
end

function C=wrapExpectedBayesianFIMConstrained(NL,NR,tauA,tauB,M,fnames,tobs_mat, ...
                                Yobs_mat,model,method,settings)
[~,tmeas_prop]=getSamplingSchedules(NL,NR,tauA,tauB);
C=expectedBayesianFIM(M,fnames,tmeas_prop,tobs_mat,Yobs_mat,model,method,settings);
end