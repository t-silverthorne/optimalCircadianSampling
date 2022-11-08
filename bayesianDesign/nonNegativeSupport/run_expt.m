clear
figure(1)
clf
figure(2)
clf
figure(3)
clf
addpath('utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

nexp=4;
settings.model='cosinorOneFreq';
settings.repar=false; % dont reparameterize in optimization step
plot_error=true;
plot_posterior=true;

rng(12)
tic
% set up time grid
NL=3;
NR=5;
tauA=0;
tauB=1/3;
settings.dT=0.2;
settings.Ntot=NL+NR;
settings.verbose=false;
settings.run_gpu=false;
settings.Tprior=20;
settings.T=50;
settings.sig1=1;
settings.sig2=settings.sig1;
settings.sig3=2*settings.sig1;
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

run_default=false;
if run_default
    nexp=4;
    settings.model='cosinorOneFreq';
    settings.repar=false; % dont reparameterize in optimization step
    plot_error=true;
    plot_posterior=true;
    
    rng(12)
    tic
    % set up time grid
    NL=3;
    NR=5;
    tauA=0;
    tauB=1/3;
    settings.dT=0.2;
    settings.Ntot=NL+NR;
    settings.verbose=false;
    settings.run_gpu=false;
    settings.Tprior=20;
    settings.T=50;
    settings.sig1=1;
    settings.sig2=settings.sig1;
    settings.sig3=2*settings.sig1;
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
    ptrue=[1.3 pi 2.5]; % true parameters
end

%%
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
[amp_est,acro_est,per_est]=convertToCircularParams(coeffvalues(bestfit),settings.model);

settings.mu1=amp_est;
settings.mu2=acro_est;
settings.mu3=1/per_est;

% error plot
if plot_error
    figure(1)
    tiledlayout(3,1,'TileSpacing','tight')
    nexttile(1)
    semilogy(1,abs(amp_est-ptrue(1))/ptrue(1),'.k','MarkerSize',15)
    hold on
    
    nexttile(2)
    semilogy(1,abs(acro_est-ptrue(2))/ptrue(2),'.k','MarkerSize',15)
    hold on
    
    nexttile(3)
    semilogy(1,abs(1/per_est-ptrue(3))/ptrue(3),'.k','MarkerSize',15)
    hold on
    drawnow
end
%%
if plot_posterior
    prior=sampleTruncatedPrior(Nsamp,settings);
    for ind=2:3 % one for unif and one for non-unif
        if ind==2
            post=samplePosteriorMCMC(Nsamp,yvec_unif,tvec_unif,settings);
        elseif ind==3
            post=samplePosteriorMCMC(Nsamp,yvec_nu,tvec_nu,settings);
        end

        figure(ind)
        tiledlayout(3,1,'TileSpacing','tight')
        nexttile
        red=[.81 .1 .26];
        blue=[.36 .54 .66];
        colours=[linspace(red(1),blue(1),nexp+1)' linspace(red(2),blue(2),nexp+1)' linspace(red(3),blue(3),nexp+1)'];
        for j=1:3
            nexttile(j)
            histogram(prior(:,j),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(1,:))
            hold on
            histogram(post(:,j),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(2,:))
            drawnow
        end
    end
end
%%

for ind=2:nexp
    % find strategy for next set of measurements
    sampstrat=particleswarm(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,tobs_nu_MAT,Yobs_nu_MAT,settings), ...
                                    2,[0,0],[1,1],opts);
    [~,mt_nu]=getSamplingSchedules(NL,NR,sampstrat(1),sampstrat(2));
    fprintf('%f\n',expectedBayesianFIM(M,mt_unif,tobs_unif_MAT,Yobs_unif_MAT,settings))
    
    % measure with new strategy
    switch settings.model 
        case 'cosinorOneFreq'
            Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
            Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));
    end    
    Yobs_unif_MAT=vertcat(Yobs_unif,Yobs_unif_MAT);
    Yobs_nu_MAT=vertcat(Yobs_nu,Yobs_nu_MAT);
    tobs_unif_MAT=vertcat(mt_unif,tobs_unif_MAT);
    tobs_nu_MAT=vertcat(mt_nu,tobs_nu_MAT);
    tvec_unif=reshape(tobs_unif_MAT,1,numel(tobs_unif_MAT));
    yvec_unif=reshape(Yobs_unif_MAT,1,numel(Yobs_unif_MAT));
    tvec_nu=reshape(tobs_nu_MAT,1,numel(tobs_nu_MAT));
    yvec_nu=reshape(Yobs_nu_MAT,1,numel(Yobs_nu_MAT));
    
    % regression with unif and non-unif measurements
    bestfit_unif=multiStartRegression(tvec_unif,yvec_unif,settings.model);
    [amp_est_unif,acro_est_unif,per_est_unif]=convertToCircularParams(coeffvalues(bestfit_unif),settings.model);
    bestfit_nu=multiStartRegression(tvec_nu,yvec_nu,settings.model);
    [amp_est_nu,acro_est_nu,per_est_nu]=convertToCircularParams(coeffvalues(bestfit_nu),settings.model);
    
    if plot_error
        figure(1)
        nexttile(1)
        semilogy(ind,abs(amp_est_unif-ptrue(1))/ptrue(1),'.k','MarkerSize',15)
        semilogy(ind,abs(amp_est_nu-ptrue(1))/ptrue(1),'.b','MarkerSize',15)
        hold on
        
        nexttile(2)
        semilogy(ind,abs(acro_est_unif-ptrue(2))/ptrue(2),'.k','MarkerSize',15)
        semilogy(ind,abs(acro_est_nu-ptrue(2))/ptrue(2),'.b','MarkerSize',15)
        hold on
        
        nexttile(3)
        semilogy(ind,abs(1/per_est_unif-ptrue(3))/ptrue(3),'.k','MarkerSize',15)
        semilogy(ind,abs(1/per_est_nu-ptrue(3))/ptrue(3),'.b','MarkerSize',15)
        hold on
    end

    if plot_posterior
        post_unif=samplePosteriorMCMC(Nsamp,yvec_unif,tvec_unif,settings);
        post_nu=samplePosteriorMCMC(Nsamp,yvec_nu,tvec_nu,settings);
        figure(2)
        for j=1:3
            nexttile(j)
            histogram(post_unif(:,j),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(ind,:))
            hold on
            drawnow
        end 
        figure(3)
        for j=1:3
            nexttile(j)
            histogram(post_nu(:,j),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(ind,:))
            hold on
            drawnow
        end 

    end
    
end
toc


%%
for ind=2:3
    figure(ind)
    for j=1:3
        nexttile(j)
        xline(ptrue(j))
    end
end
%%
% settings.repar=false;
% options=optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',1e-2,'TolX',1e-2);
% sampstrat=fminsearch(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,tobs_nu_MAT,Yobs_nu_MAT,settings),[0,1/2],options);   
% %%
% simulannealbnd(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,tobs_nu_MAT,Yobs_nu_MAT,settings),[0,1/2])
% %%
%opts=optimoptions("patternsearch",'Display','iter','PlotFcn',@psplotbestf,'UseParallel',true);
 %   patternsearch(@(x) -wrapExpectedBayesianFIM(NL,NR,x(1),x(1)+x(2),M,tobs_nu_MAT,Yobs_nu_MAT,settings),[0,1/2])
function stop=pswmyfun(optimValues,state)
stop=false;
scatter3(optimValues.swarm(:,1),optimValues.swarm(:,2),optimValues.swarmfvals)
zlim([-20,-12])
end
%%
function C=wrapExpectedBayesianFIM(NL,NR,tauA,tauB,M,tobs_mat,Yobs_mat,settings)
if settings.repar
    tauA=1+0.5*tanh(tauA);
    tauB=1+0.5*tanh(tauB);
end
[~,tmeas_prop]=getSamplingSchedules(NL,NR,tauA,tauB);
C=expectedBayesianFIM(M,tmeas_prop,tobs_mat,Yobs_mat,settings);
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
