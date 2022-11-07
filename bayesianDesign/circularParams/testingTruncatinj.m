close all
%%
tic
Nsamp=1e6;
settings.model='cosinorOneFreq';
settings.run_gpu=true;
settings.Tprior=50;
settings.sig1=1;
settings.sig2=settings.sig1;
settings.sig3=settings.sig1;
settings.mu1=1;
settings.mu2=0;
settings.mu3=2;
pr=sampleTruncatedPrior(Nsamp,settings);
toc

histogram(pr(:,1),floor(sqrt(Nsamp)),'Normalization','pdf','EdgeColor','none')
hold on
%%
xlim([0 12])

%%
% idea: use inverse-transform to generate samples from proposal function in
% hit and run, since using (6.12) acceptance probability simplifies