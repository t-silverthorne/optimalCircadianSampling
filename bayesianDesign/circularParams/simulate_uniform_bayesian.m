clear
addpath('utils/') 
addpath('FIMs/') 
addpath('models/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

rng(12)
% get estimate of parameters from uniform experiment
NL=15; % set up time grid
NR=15;
tauA=0;
tauB=1/3;
Nparam=100;
settings.dT=0.2;
settings.Ntot=NL+NR;
settings.verbose=false;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);
tobs_mat_nu=[];
Yobs_mat_nu=[];
model='cosinorOneFreq'; % get true parameter value
M=getBayesianFIMcirc(NL+NR,model); % just a function nothing evaluated yet

method='test-spt';
[theta,fnames]=samplePrior(Nparam+1,model,method,settings);
ptrue=theta(1,:);
ptrue(2)=mod(ptrue(2),2*pi); 
if settings.verbose
    ptrue
end
ptrue
theta=theta(2:end,:);

switch model % simulate measurement
    case 'cosinorOneFreq'
        Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
        %Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));
end       
Yobs_unif_MAT=Yobs_unif;
tobs_unif_MAT=mt_unif;
method='ms-prior';

% use results of multistart regression to construct prior
bestfit=multiStartRegression(mt_unif,Yobs_unif,model);
[amp_est,acro_est,per_est]=convertToCircularParams(coeffvalues(bestfit),model);
amp_est=ptrue(1);
acro_est=ptrue(2);
per_est=ptrue(3);

close all
settings.amp_est=amp_est;
settings.acro_est=acro_est;
settings.per_est=per_est;
settings.sig=1 ;

settings.speed='fast';% options: slow, fast
settings.run_gpu=false;
settings.parallel_mode='vectorize';
settings.proposal_method='fixed';
settings.FIM_expectation_method='variance';
settings.batch_size=1e4;
settings.var_cut=1e-1;

prior=samplePrior(1e4,model,method,settings);
post=samplePosteriorMCMC(1e4,fnames,mt_unif,Yobs_unif,model,method,settings);

tiledlayout(3,1,'TileSpacing','tight')
nexttile
red=[.81 .1 .26];
blue=[.36 .54 .66];
nexp=5;
colours=[linspace(red(1),blue(1),nexp+1)' linspace(red(2),blue(2),nexp+1)' linspace(red(3),blue(3),nexp+1)'];
for j=1:3
    nexttile(j)
    histogram(prior(:,j),'Normalization','pdf','EdgeColor','none','FaceColor',colours(1,:))
    hold on
    histogram(post(:,j),'Normalization','pdf','EdgeColor','none','FaceColor',colours(2,:))
end
%legend
drawnow

for i=2:nexp
    switch model % simulate measurement
        case 'cosinorOneFreq'
            Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
    end   
    Yobs_unif_MAT=[Yobs_unif; Yobs_unif_MAT];
    tobs_unif_MAT=[mt_unif; tobs_unif_MAT];
    post=samplePosteriorMCMC(1e4,fnames,tobs_unif_MAT,Yobs_unif_MAT,model,method,settings);
    for j=1:3
        nexttile(j)
        histogram(post(:,j),'Normalization','pdf','EdgeColor','none','FaceColor',colours(i+1,:))
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
xlim([-2*pi 2*pi])
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'$$-2\pi$$', '$$-\pi$$', '$$0$$', '$$\pi$$','$$2\pi$$'})
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
