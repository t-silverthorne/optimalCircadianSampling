clf
clear
addpath('../utils/')
Ntimes=31; % do not change


frq1=2;
frq2=.5;
f= [frq1 frq2 -1];
A=-[frq1 frq2  0; zeros(2,3)];
b=-[Ntimes-2; zeros(2,1)];
Aeq=[1 1 0; 0 0 1; 0 0 0];
beq=[24; Ntimes;0];
lb=[1 1 1];
nui=intlinprog(f,1:3,A,b,Aeq,beq,lb)
Ntimes=ceil(nui(1:2)'*[frq1 frq2 ]')


% A=[4 2 0.5; 1 1 1; 1 0 0];
% % b30=[Ntimes; 24; 3];
% b=[Ntimes; 24; 2];
% nui=A\b


% uniform grid 
mt_unif=linspace(0,24,Ntimes+1);
mt_unif=mt_unif(1:end-1);

% non-uniform grid
mt1=linspace(0,nui(1),frq1*nui(1)+1);
mt1=mt1(1:end-1);
mt2=linspace(nui(1),nui(1)+nui(2),ceil(frq2*nui(2))+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2];

xline(mt_nu,'-k')
xticks([0 3 8 24])

% Test condition number
% clf
% Nsamples=5;
% nreps=10000;
% zts=repmat(mt_nu,1,Nsamples);
% 
% pers=logspace(log10(.9),log10(1.1),101);
% cs=NaN(1,numel(pers));
% for i=1:numel(pers)
%     cs(i)=get_cond(zts,pers(i));
% end
% semilogy(pers,cs,'.k')
% xlim([.9 1.1])
close all

Nsamples=3;
nreps=1e3;

per1=8;
per2=4;
zts_unif = repmat(mt_unif,1,Nsamples);
zts_nu   = repmat(mt_nu,1,Nsamples);
a0=0; a1=0; a2=0; a3=.5; a4=.1;

avec=[a0 a1 a2 a3 a4 per1 per2];

sig=.1;
get_Xdat = @(zts) a0+a1*sin(2*pi*zts/per1)+a2*cos(2*pi*zts/per1) + a3*sin(2*pi*zts/per2)+a4*cos(2*pi*zts/per2) + sig*randn(nreps,numel(zts));

Xdat_unif=get_Xdat(zts_unif);
Xdat_nu=get_Xdat(zts_nu);
% plot(mt_nu,Xdat_nu(1,1:numel(mt_nu)),'.k')


per1guess=per1*(1+0*rand);
per2guess=per2*(1+0*rand);

results_unif=fit_biharmonic_cosinor(Xdat_unif,zts_unif,per1guess,per2guess);
results_nu=fit_biharmonic_cosinor(Xdat_nu,zts_nu,per1guess,per2guess);
% tiledlayout('flow')
% for i=1:size(results_unif,1)
%     nexttile
%     histogram(results_unif(i,:),ceil(sqrt(nreps)))
%     hold on
%     xline(avec(i),'LineWidth',3)
%     histogram(results_nu(i,:),ceil(sqrt(nreps)))
% end

% 
% clf
% ss=randsample(1:numel(nreps),1);
% tiledlayout(2,1)
% nexttile
% nonlinfit(zts_unif, Xdat_unif(ss,:));
% nexttile
% nonlinfit(zts_nu, Xdat_nu(ss,:));

res_nu=cell(nreps,1);
res_unif=cell(nreps,1);
tic
parfor ii=1:nreps
    res_nu{ii}=nonlinfit(zts_nu, Xdat_nu(ii,:));
    res_unif{ii}=nonlinfit(zts_unif, Xdat_unif(ii,:));
end
toc
%% Extract parameter estimates
for ii=1:numel(res_nu)
    res_nu{ii}=[res_nu{ii}.a0 res_nu{ii}.a1 res_nu{ii}.a2 res_nu{ii}.a3 res_nu{ii}.a4 res_nu{ii}.per1 res_nu{ii}.per2];
end
res_nu=cell2mat(res_nu);
for ii=1:numel(res_unif)
    res_unif{ii}=[res_unif{ii}.a0 res_unif{ii}.a1 res_unif{ii}.a2 res_unif{ii}.a3 res_unif{ii}.a4 res_unif{ii}.per1 res_unif{ii}.per2];
end
res_unif=cell2mat(res_unif);
%%
close all
tiledlayout('flow')
for i=1:size(res_unif,2)
    min_bin =min(vertcat(res_unif(:,i),res_nu(:,i)));
    max_bin =max(vertcat(res_unif(:,i),res_nu(:,i)));
    nexttile
    histogram(res_unif(:,i),ceil(sqrt(nreps)),BinLimits=[min_bin max_bin])
    hold on
    xline(avec(i),'LineWidth',3)
    histogram(res_nu(:,i),ceil(sqrt(nreps)),BinLimits=[min_bin max_bin])
end
%%
% %% Compare mean variances
% clf
% nreps=1e5;
% Vvec_unif=NaN(5,nreps);
% Vvec_nu=NaN(5,nreps);
% Nsamp=100;
% 
% 
% fprintf('\n')
% per1=24;
% per2=1;
% a0=0;
% a1=.5;
% a2=.5;
% a3=2;
% a4=.5;
% sig=.5;
% get_Xdat = @(t,Nsamp) a0+a1*sin(2*pi*t/per1)+a2*cos(2*pi*t/per1) + a3*sin(2*pi*t/per2)+a4*cos(2*pi*t/per2) + sig*randn(Nsamp,numel(t));
% 
% for ii=1:nreps
%     
%     per1guess=per1*(1+.3*rand);
%     per2guess=per2*(1+.3*rand);
%     Vvec_unif(:,ii)=var(fit_biharmonic_cosinor(get_Xdat(mt_unif,Nsamp),mt_unif,per1guess,per2guess),[],2);
%     % fprintf('Uniform: \n')
%     % for i=1:numel(V)
%     %     fprintf('%.16e \n',V(i))
%     % end
%     
%     Vvec_nu(:,ii)=var(fit_biharmonic_cosinor(get_Xdat(mt_nu,Nsamp),mt_nu,per1guess,per2guess),[],2);
%     % fprintf('Non-uniform \n')
%     % for i=1:numel(V)
%     %     fprintf('%.16e \n',V(i))
%     % end
% end
% fprintf('\n')
% fprintf('unif:     %.16e \n',mean(Vvec_unif,2))
% fprintf('non-unif: %.16e \n',mean(Vvec_nu,2))
%%

%%
% plot(mt_unif,get_Xdat(mt_unif,10),'.k')
function cc=get_cond(zts,per1,per2)

x1=sin(2*pi*zts/per1);
x2=cos(2*pi*zts/per1);
x0=ones(1,numel(zts));
X= [x0' x1' x2'];

cc=cond(X'*X);
end
