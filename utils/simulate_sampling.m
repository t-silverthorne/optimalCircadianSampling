function [res_now,gof_now]=simulate_sampling(param_true,mt_now,nreps)
% simulate sampling from biharmonic model
regularize=true; % cutoff window used in linear regression
if nargin<3
    nreps=1000;
end
warning('off','MATLAB:nearlySingularMatrix')
pctRunOnAll warning off
Nsamples=1;
beta0=param_true.beta0;
beta1=param_true.beta1;
beta2=param_true.beta2;
beta3=param_true.beta3;
beta4=param_true.beta4;
per1 =param_true.T1;
per2 =param_true.T2;
sig=.1; % noise level

zts_now   = repmat(mt_now,1,Nsamples);

% for generating data
get_Xdat = @(zts) beta0+beta1*sin(2*pi*zts/per1)+beta2*cos(2*pi*zts/per1) + beta3*sin(2*pi*zts/per2)+beta4*cos(2*pi*zts/per2) + sig*randn(nreps,numel(zts));
Xdat_now=get_Xdat(zts_now); % sample on non-uniform grid

res_now=cell(nreps,1);
gof_now=cell(nreps,1);

Xdat_now=parallel.pool.Constant(Xdat_now);
if nreps>1
    parfor ii=1:nreps
        [res_now{ii},gof_now{ii}]=nonlinfit_grid_fast(zts_now, Xdat_now.Value(ii,:),regularize);
    end
else
    [res_now{1},gof_now{1}]=nonlinfit_grid_fast(zts_now, Xdat_now.Value(1,:),regularize);
end
warning('on','MATLAB:nearlySingularMatrix')
pctRunOnAll warning on
end