%% For comparing performance of code with GPU acceleration

clear
addpath('utils/') 
addpath('FIMs/') 
addpath('models/')

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
        Yobs_nu=cosinorOneFreq(mt_nu,getTheta(ptrue,fnames))+randn(1,numel(mt_nu));
end       
tobs_mat_nu=mt_nu;
Yobs_mat_nu=Yobs_nu;

% use results of multistart regression to construct prior
bestfit=multiStartRegression(mt_unif,Yobs_unif,model);
[amp_est,acro_est,per_est]=convertToCircularParams(coeffvalues(bestfit),model);
settings.amp_est=amp_est;
settings.acro_est=acro_est;
settings.per_est=per_est;

%%
Nbench=1e4;
settings.parallel_mode='vectorize';
method='test-spt'; % TODO: make this Ms prior
settings.proposal_method='fixed'; % options: fixed, iterative

settings.speed='slow';% options: slow, fast
settings.run_gpu=false;
tic
samplePosteriorMCMC(Nbench,fnames,tobs_mat_nu, ...
                Yobs_mat_nu,model,method,settings);
timenow=toc;
fprintf('\nRecursive      : %f\n',timenow)


settings.speed='slow';% options: slow, fast
settings.run_gpu=true;
tic
samplePosteriorMCMC(Nbench,fnames,tobs_mat_nu, ...
                Yobs_mat_nu,model,method,settings);
timenow=toc;
fprintf('Recursive w GPU: %f\n',timenow)


settings.speed='fast';% options: slow, fast
settings.run_gpu=false;
tic
samplePosteriorMCMC(Nbench,fnames,tobs_mat_nu, ...
                Yobs_mat_nu,model,method,settings);
timenow=toc;
fprintf('Vectorize      : %f\n',timenow)


settings.speed='fast';% options: slow, fast
settings.run_gpu=true;
tic
samplePosteriorMCMC(Nbench,fnames,tobs_mat_nu, ...
                Yobs_mat_nu,model,method,settings);
timenow=toc;
fprintf('Vectorize w GPU: %f\n',timenow)


