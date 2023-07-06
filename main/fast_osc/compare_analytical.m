addpath('../../guelph_talk/utils_core/')
clear all
Nmeas=4;
mt=linspace(0,1,Nmeas+1);
mt=mt(1:end-1);
nMC=50;
Nres=1e4;

Amp=9;
freq=5.4;
maxit=1e3;
acro=0;
p.Amp=Amp;
p.freq=freq;
p.acro=acro;
X=constructX(mt,p);


% verify chisq(lambda)
csq    = cos(2*pi*freq*mt-acro);
lambdapap=Amp^2*csq*csq';

mu=[0; Amp*cos(acro); Amp*sin(acro)];
lambda=mu'*((X'*X)\mu);
%
Y=Amp*cos(2*pi*freq*mt-acro)+randn(Nres,Nmeas);

betas_obs = pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y)));

close all
num=sum(betas_obs.*((X'*X)\betas_obs),1)';


x=linspace(min(num),max(num),1e3);
histogram(num,'Normalization','probability')
hold on
plot(x,ncx2pdf(x,3,lambda),'-k')

%%
% pwr_exact=getPowerExact(Amp,acro,freq,mt);
% pwr_est=0;
% for kk=1:nMC
% 
% Y=Amp*cos(2*pi*freq*mt-acro)+randn(Nres,Nmeas);
% 
% X=constructX(mt,p);
% betas_obs = pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y)));
% fits_obs  = pagetranspose(pagemtimes(X,betas_obs));
% SSres_obs = sqrt(sum((fits_obs-Y).^2,2));
% 
% YI=NaN(size(Y,1),size(Y,2),factorial(Nmeas));
% all_perms=perms(1:Nmeas);
% for ii=1:factorial(Nmeas)
%     YI(:,:,ii)=Y(:,all_perms(ii,:));
% end
% betas     = pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
% fits      = pagetranspose(pagemtimes(X,betas));
% SSres     = sqrt(sum((fits-YI).^2,2));
% 
% 
% bin_mat   = SSres<SSres_obs;
% pwr_est=pwr_est+sum(sum(bin_mat,3)/factorial(Nmeas) <.05)/Nres;
% end
% pwr_est/nMC