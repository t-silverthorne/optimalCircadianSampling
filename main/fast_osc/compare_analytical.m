addpath('../../guelph_talk/utils_core/')
clear all
Nmeas=4;
mt=linspace(0,1,Nmeas+1);
mt=mt(1:end-1);
nMC=2;
Nres=5e4;

Amp=1;
freq=3;
acro=0;
p.Amp=Amp;
p.freq=freq;
p.acro=acro;
X=constructX(mt,p);


% verify chisq(lambda)
csq       = cos(2*pi*freq*mt-acro);
lambdapap = Amp^2*csq*csq';


%% simulate

Y=Amp*cos(2*pi*freq*mt-acro)+randn(Nres,Nmeas);

ft=fit_cosinor_model(Y,mt,1/freq)

sum(ft.pval<.05)/Nres

%%

% mu=[0; Amp*sin(acro); Amp*cos(acro)];
% lambda=mu'*((X'*X)*mu);


% Y=Amp*cos(2*pi*freq*mt-acro)+randn(Nres,Nmeas);

% betas_obs = pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y)));
% 
% close all
% num=sum(betas_obs.*((X'*X)*betas_obs),1)';
% 
% x=linspace(min(num),max(num),1e3);
% tiledlayout(2,1)
% % histogram(num,'Normalization','pdf')
% % hold on
% % plot(x,ncx2pdf(x,2,lambda),'-k', 'linewidth',2)
% 
% % verify F distribution
% fits_obs  = pagetranspose(pagemtimes(X,betas_obs));
% SSres_obs = sum((fits_obs-Y).^2,2);
% histogram((Nmeas-3)*num./SSres_obs/2,0:2:1e3,'Normalization','pdf')
% x=linspace(0,600,1e3);
% hold on
% plot(x,ncfpdf(x,2,Nmeas-3,lambda),'-k', 'linewidth',2)

pwr_exact=getPowerExact(Amp,acro,freq,mt);
pwr_est=0;

for nn=1:nMC
Y=Amp*cos(2*pi*freq*mt-acro)+randn(Nres,Nmeas);

X=constructX(mt,p);
betas_obs = pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y)));
fits_obs  = pagetranspose(pagemtimes(X,betas_obs));
SSres_obs  = sum((fits_obs-Y).^2,2);

X0 = ones(Nmeas,1);
betas_obs_0 = pagemldivide(X0'*X0,pagemtimes(X0',pagetranspose(Y)));
fits_obs_0  = pagetranspose(pagemtimes(X0,betas_obs_0));
SSres_obs0  = sum((fits_obs_0-Y).^2,2);
Fobs        = (SSres_obs0-SSres_obs)./SSres_obs;


YI=NaN(size(Y,1),size(Y,2),factorial(Nmeas));
all_perms=perms(1:Nmeas);
if Nmeas<9
    for ii=1:factorial(Nmeas)
        YI(:,:,ii)=Y(:,all_perms(ii,:));
        Nperm=factorial(Nmeas);
    end
end
betas     = pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits      = pagetranspose(pagemtimes(X,betas));
SSres     = sum((fits-YI).^2,2);

betas_0 = pagemldivide(X0'*X0,pagemtimes(X0',pagetranspose(YI)));
fits_0  = pagetranspose(pagemtimes(X0,betas_0));
SSres_0 = sum((fits_0-Y).^2,2);
FI      = (SSres_0-SSres)./SSres;

%bin_mat   = SSres<SSres_obs;
bin_mat  = Fobs < FI;
pwr_est=pwr_est+sum(sum(bin_mat,3)/Nperm <.05)/Nres;
end
pwr_exact
pwr_est/nMC