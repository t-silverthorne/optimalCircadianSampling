%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks that exact and approximate methods agree
addpath('../utils')

% parameters
Nmeas=4; % strongest when number of measurements is small
mt=linspace(0,1,Nmeas+1);
mt=mt(1:end-1);
Nres=5e4;
Amp=2+rand;
freq=2+rand;
acro=rand*2*pi;
p.Amp=Amp;
p.freq=freq;
p.acro=acro;
X=constructX(mt,p);
fprintf('\nParameters: \n')
fprintf('   Nmeas = %d \n',Nmeas)
fprintf('   Nres  = %d \n',Nres)
fprintf('   Amp   = %.4f \n',Amp)
fprintf('   freq  = %.4f \n',freq)
fprintf('   acro  = %.4f \n\n',acro)

% get exact power
pwr_exact=getPowerExact(Amp,acro,freq,mt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo estimates (F-test and error test)
Y=Amp*cos(2*pi*freq*mt-acro)+randn(Nres,Nmeas);

X=constructX(mt,p);
betas_obs   = pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y)));
fits_obs    = pagetranspose(pagemtimes(X,betas_obs));
SSres_obs   = sum((fits_obs-Y).^2,2);

YI=NaN(size(Y,1),size(Y,2),factorial(Nmeas));
all_perms=perms(1:Nmeas);
if Nmeas<9
    for ii=1:factorial(Nmeas)
        YI(:,:,ii)=Y(:,all_perms(ii,:));
        Nperm=factorial(Nmeas);
    end
else
    pwr_est=NaN;
end

betas     = pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits      = pagetranspose(pagemtimes(X,betas));
SSres     = sum((fits-YI).^2,2);

bin_mat  = SSres_obs < SSres;
pwr_est_T =sum(sum(bin_mat,3)/Nperm <.05)/Nres;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test using cosinor model
ft=fit_cosinor_model(Y,mt,1/freq);
pwr_est_f=mean(ft.pval<.05);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print results
fprintf('Exact power estimate:   %.4f \n',pwr_exact)
fprintf('MC power estimate:      %.4f \n',pwr_est_T)
fprintf('F-test power estimate:  %.4f \n\n',pwr_est_f)
