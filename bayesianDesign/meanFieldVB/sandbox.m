clear
figure(1)
clf
figure(2)
clf
figure(3)
clf
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

addpath('../nonNegativeSupport/utils/')

nexp=4;
settings.model='cosinorOneFreq';
settings.repar=false; % dont reparameterize in optimization step
plot_error=true;
plot_posterior=true;

rng(12)
%% set up time grid and monte carlo integration
NL=3;
NR=5;
tauA=0;
tauB=1/2;
settings.dT=0.2;
settings.Ntot=NL+NR;
settings.verbose=false;
settings.run_gpu=false;
settings.Tprior=20;
settings.T=50;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

% for terminating Monte Carlo integration of FIM
settings.var_cut=1e-1;
settings.batch_size=1e4;
settings.FIM_expectation_method='variance';
Nsamp=8*settings.batch_size;

M=getBayesianFIMcirc(NL+NR,settings.model); % just a function nothing evaluated yet

opts = optimoptions(@particleswarm, ... %'HybridFcn',@fminsearch,...
														 'Display','iter', ...    
														 'SwarmSize',16, ...
														 'MaxIter',10, ... % TODO make realistic
														 'InitialSwarmSpan',1/2, ...
														 'FunctionTolerance',1e-2, ... % might be too high
														 'UseParallel',true);%,'PlotFcn',@pswmyfun); % TODO make true


% simulate measurement both start with uniform data

ptrue=[2.7 pi 4.5]; % true parameters
ptruelin=[ptrue(1)*cos(ptrue(2)),ptrue(1)*sin(ptrue(2)),ptrue(3)]; % linear form of true

fnames={'A1','phi1','f1'};
switch settings.model % simulate measurement
    case 'cosinorOneFreq'
        Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
        %Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));
end       
Yobs_unif_MAT=Yobs_unif;
tobs_unif_MAT=mt_unif;
Yobs_nu_MAT=Yobs_unif;
tobs_nu_MAT=mt_unif;

% reshape measurements for evaluating posterior distribution
tvec_unif=reshape(tobs_unif_MAT,1,numel(tobs_unif_MAT));
yvec_unif=reshape(Yobs_unif_MAT,1,numel(Yobs_unif_MAT));
tvec_nu=reshape(tobs_nu_MAT,1,numel(tobs_nu_MAT));
yvec_nu=reshape(Yobs_nu_MAT,1,numel(Yobs_nu_MAT));

% multistart regression to get best fit
bestfit=multiStartRegression(mt_unif,Yobs_unif,settings.model);

% use bestfit to get 
settings.muA=bestfit.a2;
settings.sigA=5;
settings.muB=bestfit.a1; % coeff on sin term
settings.sigB=5;
settings.muf=1/bestfit.per1;
settings.sigf=5;

settings.Lcut=25;

%% variational Bayes
tic

[qA_nu,qB_nu,qf_nu,pars_nu,qfdom_nu,Zf0_nu,Zf_nu]=updateVariationalBayes(tvec_nu,yvec_nu,settings);
[qA_unif,qB_unif,qf_unif,pars_unif,qfdom_unif,Zf0_unif,Zf_unif]=updateVariationalBayes(tvec_unif,yvec_unif,settings);
toc

figure(2)
fvals=0:.01:settings.Lcut;
Avals=-10:.01:10;
Bvals=-10:.01:10;

tiledlayout(3,1)
nexttile(1)
plot(Avals,qA_unif(Avals))
hold on
xline(ptruelin(1))

nexttile(2)
plot(Bvals,qB_unif(Bvals))
hold on
xline(ptruelin(2))

nexttile(3)
plot(fvals,qf_unif(fvals))
hold on
xline(ptruelin(3))

figure(3)
tiledlayout(3,1)
nexttile(1)
plot(Avals,qA_nu(Avals))
hold on
xline(ptruelin(1))

nexttile(2)
plot(Bvals,qB_nu(Bvals))
hold on
xline(ptruelin(2))

nexttile(3)
plot(fvals,qf_nu(fvals))
hold on
xline(ptruelin(3))


%% sample posterior distributions

% sample Gaussians to get linear params
% samp_nu=randn(Nsamp,2);
% samp_nu(:,1)=samp_nu(:,1)*sigA_nu+muA_nu;
% samp_nu(:,2)=samp_nu(:,2)*sigB_nu+muB_nu;

% use importance sampling to get period params











