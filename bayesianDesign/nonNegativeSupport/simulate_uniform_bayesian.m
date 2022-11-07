clear
addpath('utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

rng(12)
% setup
NL=5; % set up time grid
NR=10;
tauA=0;
tauB=1/3;
Nparam=100;
Nsamp=1e4

settings.dT=0.2;
settings.Ntot=NL+NR;
settings.verbose=false;
settings.run_gpu=false;
settings.Tprior=50;
settings.T=300;
settings.sig1=1;
settings.sig2=settings.sig1;
settings.sig3=10*settings.sig1;

[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

tobs_mat_nu=[];
Yobs_mat_nu=[];
settings.model='cosinorOneFreq'; % get true parameter value
M=getBayesianFIMcirc(NL+NR,settings.model); % just a function nothing evaluated yet
% simulate measurement
ptrue=[10 pi 2];
fnames={'A1','phi1','f1'};
switch settings.model % simulate measurement
    case 'cosinorOneFreq'
        Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
end       
Yobs_unif_MAT=Yobs_unif;
tobs_unif_MAT=mt_unif;
tvec=reshape(tobs_unif_MAT,1,numel(tobs_unif_MAT));
yvec=reshape(Yobs_unif_MAT,1,numel(Yobs_unif_MAT));

% get prior
% use results of multistart regression to construct prior
bestfit=multiStartRegression(mt_unif,Yobs_unif,settings.model);
[amp_est,acro_est,per_est]=convertToCircularParams(coeffvalues(bestfit),settings.model);
amp_est
acro_est
per_est
f1_est=1/per_est;
settings.mu1=amp_est;
settings.mu2=acro_est;
settings.mu3=1/per_est;

prior=sampleTruncatedPrior(Nsamp,settings);
histogram(prior(:,3),80,'EdgeColor','none')

post=samplePosteriorMCMC(Nsamp,yvec,tvec,settings);

tiledlayout(3,1,'TileSpacing','tight')
nexttile
red=[.81 .1 .26];
blue=[.36 .54 .66];
nexp=5;
colours=[linspace(red(1),blue(1),nexp+1)' linspace(red(2),blue(2),nexp+1)' linspace(red(3),blue(3),nexp+1)'];
for j=1:3
    nexttile(j)
    histogram(prior(:,j),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(1,:))
    hold on
    histogram(post(:,j),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(2,:))
end
%legend
drawnow

for i=2:nexp
    switch settings.model % simulate measurement
        case 'cosinorOneFreq'
            Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
    end   
    Yobs_unif_MAT=[Yobs_unif; Yobs_unif_MAT];
    tobs_unif_MAT=[mt_unif; tobs_unif_MAT];
    post=samplePosteriorMCMC(Nsamp,yvec,tvec,settings);
    for j=1:3
        nexttile(j)
        histogram(post(:,j),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(i+1,:))
        hold on
    end
    drawnow
end
for j=1:3
    nexttile(j)
    xline(ptrue(j))
end
hold off

nexttile(1)
xlabel('amp')

nexttile(2)
xlim([0 2*pi])
xticks([0 pi 2*pi])
xticklabels({'$$0$$', '$$\pi$$','$$2\pi$$'})
xlabel('acro')

nexttile(3)
xlabel('period')

nvals=0:1:(nexp+1);
h=colorbar;
caxis([0 nexp+1])
h.Limits=[0 nexp+1];

h.Ticks=nvals+.5;

set(h,'TickLabelInterpreter','latex')
colormap(colours)
h.Layout.Tile='east'

h.TickLabels=num2cell(nvals+1);

h.Label.Interpreter='latex'
h.Label.String='experiment';
