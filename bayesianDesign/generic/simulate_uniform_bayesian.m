addpath('utils/') 
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nhist=1e5

rng(12)
% get estimate of parameters from uniform experiment
NL=4; % set up time grid
NR=8;
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
settings.nondeg=true;
fnames={'A1','phi1','f1'};
method='test-spt';

theta=[10 0 1/0.6]; % TODO make prior generation better
ptrue=theta(1,:);
ptrue(2)=mod(ptrue(2),2*pi); 
if settings.verbose
    ptrue
end
ptrue(1)=sqrt(ptrue(1));
makeDeg(makeNonDeg(ptrue,model),model)
%%

switch model % simulate measurement
    case 'cosinorOneFreq'
        if settings.nondeg
            Yobs_unif=cosinorOneFreq(mt_unif,getTheta(makeNonDeg(ptrue,model),fnames))+randn(1,numel(mt_unif));
        else
            Yobs_unif=cosinorOneFreq(mt_unif,getTheta(ptrue,fnames))+randn(1,numel(mt_unif));
        end
        
end       

Yobs_unif_MAT=Yobs_unif;
tobs_unif_MAT=mt_unif;



% use results of multistart regression to construct prior
bestfit=multiStartRegression(mt_unif,Yobs_unif,model);
[amp_est,acro_est,per_est]=convertToCircularParams(coeffvalues(bestfit),model);
freq_est=1/per_est;
deg_pvec=makeDeg([amp_est acro_est freq_est],model);


close all
settings.amp_est=deg_pvec(1);
settings.acro_est=deg_pvec(2);
settings.freq_est=deg_pvec(3);
settings.sig=2 ;

settings.run_gpu=true;
settings.proposal_method='fixed';


method='ms-prior';
prior=samplePrior(Nhist,model,method,settings);
prior=real(cell2mat(arrayfun(@(ii)  makeNonDeg(prior(ii,:),model)',1:size(prior,1),'UniformOutput',false))');

post=samplePosteriorMCMC(Nhist,fnames,mt_unif,Yobs_unif,model,method,settings);
post=real(cell2mat(arrayfun(@(ii)  makeNonDeg(post(ii,:),model)',1:size(post,1),'UniformOutput',false))');

tiledlayout(3,1,'TileSpacing','tight')
nexttile
red=[.81 .1 .26];
blue=[.36 .54 .66];
nexp=5;
colours=[linspace(red(1),blue(1),nexp+1)' linspace(red(2),blue(2),nexp+1)' linspace(red(3),blue(3),nexp+1)'];
for j=1:3
    nexttile(j)
    histogram(prior(:,j),floor(sqrt(Nhist)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(1,:),'FaceAlpha',0.6)
    hold on
    histogram(post(:,j),floor(sqrt(Nhist)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(2,:),'FaceAlpha',0.1)
end
%legend
drawnow


for i=2:nexp
    switch model % simulate measurement
        case 'cosinorOneFreq'
            Yobs_unif=cosinorOneFreq(mt_unif,getTheta(makeNonDeg(ptrue,model),fnames))+randn(1,numel(mt_unif));
    end   
    Yobs_unif_MAT=[Yobs_unif; Yobs_unif_MAT];
    tobs_unif_MAT=[mt_unif; tobs_unif_MAT];
    post=samplePosteriorMCMC(Nhist,fnames,tobs_unif_MAT,Yobs_unif_MAT,model,method,settings);
    post=real(cell2mat(arrayfun(@(ii)  makeNonDeg(post(ii,:),model)',1:size(post,1),'UniformOutput',false))');

    for j=1:3
        nexttile(j)
        if i==nexp
            histogram(post(:,j),floor(sqrt(Nhist)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(i+1,:),'FaceAlpha',0.6)
        else
            histogram(post(:,j),floor(sqrt(Nhist)),'Normalization','pdf','EdgeColor','none','FaceColor',colours(i+1,:),'FaceAlpha',0.1)
        end
        hold on
    end
    drawnow
end

ptruerepar=makeNonDeg(ptrue,model);
for j=1:3
    nexttile(j)
    xline(ptruerepar(j))
end
hold off
%%
nexttile(1)
xlabel('amp')

nexttile(2)
xlim([0 2*pi])
xticks([ 0 pi 2*pi])
xticklabels({ '$$0$$', '$$\pi$$','$$2\pi$$'})
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