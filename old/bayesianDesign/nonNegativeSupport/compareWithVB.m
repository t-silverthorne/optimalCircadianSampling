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
settings.Tprior=5e3;
settings.T=5e3;
settings.sig1=1;
settings.sig2=settings.sig1;
settings.sig3=2*settings.sig1;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

% for terminating Monte Carlo integration of FIM
settings.var_cut=1e-1;
settings.batch_size=1e4;
settings.FIM_expectation_method='variance';
Nsamp=8*settings.batch_size;

% simulate measurement both start with uniform data
ptrue=[2.7 pi 4.5]; % true parameters
ptruelin=[ptrue(1)*cos(ptrue(2)),ptrue(1)*sin(ptrue(2)),ptrue(3)]; % linear form of ptrue



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

settings.mu1=bestfit.a2;
settings.sig1=1;
settings.mu2=bestfit.a1; % coeff on sin term
settings.sig2=1;
settings.mu3=1/bestfit.per1;
settings.sig3=2;
%%
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
    prior=sampleLinTruncatedPrior(Nsamp,settings);
    for ind=2:3 % one for unif and one for non-unif
        if ind==2
            post=sampleLinPosteriorMCMC(Nsamp,yvec_unif,tvec_unif,settings);
        elseif ind==3
            post=sampleLinPosteriorMCMC(Nsamp,yvec_nu,tvec_nu,settings);
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



