close all
clear
settings.run_gpu=false;
Nsamp=1e7;
Nmeas=10;
mu_true=10;
settings.sigma=12;
settings.sigma0=30;
settings.mu0=1;

yvec=mu_true+randn([1,Nmeas])*settings.sigma; % data from true distribution


s0=settings.sigma0;
sigma=settings.sigma;
mu0=settings.mu0;
n=length(yvec);
ytot=sum(yvec);
mupost=(mu0/s0^2 + ytot/sigma^2)/(1/s0^2+n/sigma^2);
sigmapost=1/(1/s0^2+n/sigma^2);



model='gaussian';
S=samplePosteriorMCMC(1e3,[],[],yvec,'gaussian',[],settings);
histogram(S,'normalization','pdf')
xv=-30:.1:30;
hold on
plot(xv,normpdf(xv,mupost,sqrt(sigmapost)))