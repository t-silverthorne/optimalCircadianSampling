clf
clear
addpath('../utils/')

% example: uniform is not optimal
% rng(8)
% param.beta0=0;
% param.beta1=.1;
% param.beta2=1;
% param.beta3=1;
% param.beta4=.1;
% param.T1=3;
% param.T2=15;
% testing=true;
% n0=1;
% niter=100;

rng(5)
n0=12;
niter=200;
testing=false;
% true parameters 
param_true.beta0=0;
param_true.beta1=0;
param_true.beta2=.5;
param_true.beta3=.25;
param_true.beta4=.25;
param_true.T1=12;
param_true.T2=4;
theta_unif=linspace(0,24,10);
theta_unif=theta_unif(1:end-1);
% estimate from uniform
[res,gof]=simulate_sampling(param_true,theta_unif,1);
res{1}
%%
n_batches=5
for jj=1:n_batches
    param_est.beta0=res{1}.a0;
    param_est.beta1=res{1}.a1;
    param_est.beta2=res{1}.a2;
    param_est.beta3=res{1}.a3;
    param_est.beta4=res{1}.a4;
    param_est.T1=res{1}.per1;
    param_est.T2=res{1}.per2;
    theta_vec=0;
    while numel(theta_vec)<7
        tv_list={};
        Psi_list={};
        for ii=1:3
            [tv_list{ii},Psi_list{ii}]=run_atwood_nonlinear(param_est,niter,n0,testing);
        end
        [~,ind]=max(cell2mat(Psi_list));
        theta_vec=tv_list{ind};
    end
    theta_vec=24*sort(theta_vec)';
    [res,~]=simulate_sampling(param_est,theta_vec,1);
    res{1}
end


%%
theta_unif=0:2:24;
theta_unif=theta_unif(1:end-1);
theta_vec=run_atwood_nonlinear(param,niter,n0,testing);
numel(theta_vec);
[res,gof]=simulate_sampling(param,theta_unif,1)
theta_vec=24*theta_vec';
[res,gof]=simulate_sampling(param,theta_vec,1)
res{1}
%%
mean()

