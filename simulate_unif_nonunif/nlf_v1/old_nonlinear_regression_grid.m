clf
clear
addpath('../../utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
rng(1)

construct non-uniform grid
Ntimes=31; 
frq1=2;  % samples per hour in first region
frq2=.5; % sph in second region
f= [frq1 frq2 -1];
A=-[frq1 frq2  0; zeros(2,3)];
b=-[Ntimes-2; zeros(2,1)];
Aeq=[1 1 0; 0 0 1; 0 0 0];
beq=[24; Ntimes;0];
lb=[1 1 1];
nui=intlinprog(f,1:3,A,b,Aeq,beq,lb);
Ntimes=ceil(nui(1:2)'*[frq1 frq2 ]'); % number of samples consistent with restrictions

mt1=linspace(0,nui(1),frq1*nui(1)+1);
mt1=mt1(1:end-1);
mt2=linspace(nui(1),nui(1)+nui(2),ceil(frq2*nui(2))+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2];

% construct uniform grid 
mt_unif=linspace(0,24,Ntimes+1);
mt_unif=mt_unif(1:end-1);

% sampling options
Nsamples=3; % number of samples at each time point
nreps=1e3;  % number of times to repeat experiment

per1=8;
per2=4;
zts_unif = repmat(mt_unif,1,Nsamples);
zts_nu   = repmat(mt_nu,1,Nsamples);
a0=0; a1=.25; a2=.25; a3=.5; a4=.5;
sig=.1; % noise level

avec=[a0 a1 a2 a3 a4 per1 per2]; % for plotting

% for generating data
get_Xdat = @(zts) a0+a1*sin(2*pi*zts/per1)+a2*cos(2*pi*zts/per1) + a3*sin(2*pi*zts/per2)+a4*cos(2*pi*zts/per2) + sig*randn(nreps,numel(zts));

Xdat_unif=get_Xdat(zts_unif); % sample on uniform grid
Xdat_nu=get_Xdat(zts_nu); % sample on non-uniform grid

res_nu=cell(nreps,1);
res_unif=cell(nreps,1);
gof_nu=cell(nreps,1);
gof_unif=cell(nreps,1);


%%
% tic
% parfor ii=1:nreps
%     res_nu{ii}=nonlinfit_grid(zts_nu, Xdat_nu(ii,:));
%     res_unif{ii}=nonlinfit_grid(zts_unif, Xdat_unif(ii,:));
% end
% toc

Xdat_unif=parallel.pool.Constant(Xdat_unif);
Xdat_nu=parallel.pool.Constant(Xdat_nu);


tic
parfor ii=1:nreps
    [res_nu{ii},gof_nu{ii}]=nonlinfit_grid(zts_nu, Xdat_nu.Value(ii,:));
    [res_unif{ii},gof_unif{ii}]=nonlinfit_grid(zts_unif, Xdat_unif.Value(ii,:));
end
toc

%%
save(strcat('Nsamp_',...
num2str(Nsamples),... 
'_nreps_',...
num2str(nreps),...
'_per1_',...
num2str(per1),...
'_per2_',...
num2str(per2)))
%%
% Extract parameter estimates
for ii=1:numel(res_nu)
    res_nu{ii}=[res_nu{ii}.a0 res_nu{ii}.a1 res_nu{ii}.a2 res_nu{ii}.a3 res_nu{ii}.a4 res_nu{ii}.per1 res_nu{ii}.per2];
end
res_nu=cell2mat(res_nu);
for ii=1:numel(res_unif)
    res_unif{ii}=[res_unif{ii}.a0 res_unif{ii}.a1 res_unif{ii}.a2 res_unif{ii}.a3 res_unif{ii}.a4 res_unif{ii}.per1 res_unif{ii}.per2];
end
res_unif=cell2mat(res_unif);
