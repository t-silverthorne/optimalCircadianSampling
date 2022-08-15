addpath('../utils/')
a0=0; a1=0; a2=.5; a3=.25; a4=.25;
per1=12;
per2=4;
% a0=rand; a1=rand; a2=rand; a3=rand; a4=rand;
% per1=rand*1;
% per2=rand*10;
param.beta0=a0;
param.beta1=a1;
param.beta2=a2;
param.beta3=a3;
param.beta4=a4;
param.T1=per1/24;
param.T2=per2/24;

Nleft=8;
Nright=2;
tic
method='nonlin'% nonlin or lin
fprintf(strcat('\n ',method,'\n'))
for outer_index=1:10

    Nleft=Nleft+2;
    Nright=Nright+2;
    Ntimes=Nleft+Nright;
    
    mt1=linspace(0,8,Nleft+1);
    mt1=mt1(1:end-1);
    mt2=linspace(8,24,Nright+1);
    mt2=mt2(1:end-1);
    mt_nu=[mt1 mt2];
    
    % construct uniform grid 
    mt_unif=linspace(0,24,Ntimes+1);
    mt_unif=mt_unif(1:end-1);
    
    mt_nu=mt_nu/24;
    mt_unif=mt_unif/24;
    switch method
        case 'nonlin'
            Cnu=log(det(get_FIM_biharmonic_nonlin(mt_nu',ones(numel(mt_nu),1)/numel(mt_nu),param)));
            Cunif=log(det(get_FIM_biharmonic_nonlin(mt_unif',ones(numel(mt_unif),1)/numel(mt_unif),param)));
        case 'lin'
            Cnu=log(det(get_info_matrix_k(mt_nu',ones(numel(mt_nu),1)/numel(mt_nu),2)));
            Cunif=log(det(get_info_matrix_k(mt_unif',ones(numel(mt_unif),1)/numel(mt_unif),2)));
    end

    fprintf('NL: %d  NR: %d  Psival unif: %f Psival nu %f %d %d\n',Nleft,Nright,Cunif,Cnu,numel(mt_nu),numel(mt_unif))
end
toc
%%

%example: uniform is not optimal
rng(8)
param.beta0=0;
param.beta1=.1;
param.beta2=1;
param.beta3=1;
param.beta4=.1;
param.T1=3;
param.T2=15;
testing=true;
n0=1;
niter=50;

run_atwood_nonlinear(param,niter,n0,true)