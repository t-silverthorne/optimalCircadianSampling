clf
clear
addpath('../utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


Nleft=10;
Nright=4;

Ntimes=Nleft+Nright;

mt1=linspace(0,8,Nleft+1);
mt1=mt1(1:end-1);
mt2=linspace(8,24,Nright+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2];

% construct uniform grid 
mt_unif=linspace(0,24,Ntimes+1);
mt_unif=mt_unif(1:end-1);

% sampling options
Nsamples=1; % number of samples at each time point
nreps=1e3;  % number of times to repeat experiment

per1=12;
per2=4;
zts_unif = repmat(mt_unif,1,Nsamples);
zts_nu   = repmat(mt_nu,1,Nsamples);
a0=0; a1=0; a2=.5; a3=.25; a4=.25;
sig=.1; % noise level

avec=[a0 a1 a2 a3 a4 per1 per2]; % for plotting

% for generating data
get_Xdat = @(zts) a0+a1*sin(2*pi*zts/per1)+a2*cos(2*pi*zts/per1) + a3*sin(2*pi*zts/per2)+a4*cos(2*pi*zts/per2) + sig*randn(nreps,numel(zts));

Xdat_unif=get_Xdat(zts_unif); % sample on uniform grid
Xdat_nu=get_Xdat(zts_nu); % sample on non-uniform grid

%%
tic
nonlinfit(zts_nu, Xdat_nu(1,:), 9,20)
toc
%%
tic
nonlinfit_grid(zts_nu, Xdat_nu(1,:))
toc