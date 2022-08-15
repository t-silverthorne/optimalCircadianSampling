clf
clear
addpath('../utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% non-uniform grid
Ntimes=31; % Ideal number of samples
frq1=2; % samples per hour in first region
frq2=.5;
f= [frq1 frq2 -1];
A=-[frq1 frq2  0; zeros(2,3)];
b=-[Ntimes-2; zeros(2,1)];
Aeq=[1 1 0; 0 0 1; 0 0 0];
beq=[24; Ntimes;0];
lb=[1 1 1];
nui=intlinprog(f,1:3,A,b,Aeq,beq,lb);
Ntimes=ceil(nui(1:2)'*[frq1 frq2 ]'); % number of samples consistent with restrictions

% finish non-uniform grid
mt1=linspace(0,nui(1),frq1*nui(1)+1);
mt1=mt1(1:end-1);
mt2=linspace(nui(1),nui(1)+nui(2),ceil(frq2*nui(2))+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2];

% uniform grid 
mt_unif=linspace(0,24,Ntimes+1);
mt_unif=mt_unif(1:end-1);

% show grids
% xline(mt_nu,'-k')
% xticks([0 3 8 24])
% close all

% Number of samples at each time points, number of times to repeat
% experiment
Nsamples=3;
nreps=1e3;

% per1=8;
% per2=4;
% zts_unif = repmat(mt_unif,1,Nsamples);
% zts_nu   = repmat(mt_nu,1,Nsamples);
% a0=0; a1=0; a2=0; a3=.5; a4=.1;

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

per1guess=per1*(1+0*rand); % don't use this
per2guess=per2*(1+0*rand);

results_unif=fit_biharmonic_cosinor(Xdat_unif,zts_unif,per1guess,per2guess);
results_nu=fit_biharmonic_cosinor(Xdat_nu,zts_nu,per1guess,per2guess);
%%


% show histograms for linear regression problem
tiledlayout(1,10,'Padding','tight','TileSpacing','tight')
for i=1:size(results_unif,1)
    min_bin =min(horzcat(results_unif(i,:),results_nu(i,:)));
    max_bin =max(horzcat(results_unif(i,:),results_nu(i,:)));
    nexttile([1 2])
    histogram(results_unif(i,:),ceil(sqrt(nreps)),'LineStyle','none',BinLimits=[min_bin max_bin])
    hold on
    histogram(results_nu(i,:),ceil(sqrt(nreps)),'LineStyle','none',BinLimits=[min_bin max_bin])
    xline(avec(i),'-k','LineWidth',1,'Alpha',1,'Color','black')
    ylim([0 120])
end

nexttile(1)
ylabel('count','interpreter','latex')
title('$\beta_0$','Interpreter','latex')

for i=[3 5 7 9]
    nexttile(i)
    set(gca,'YTickLabel',[]);
    %xlabel('count','interpreter','latex')
end

nexttile(5)
xlabel('linear regression estimate','interpreter','latex')
legend({'uniform','non-uniform','true value'},'Location','southoutside','NumColumns',3)

nexttile(3)
title('$\beta_1$','Interpreter','latex')
nexttile(5)
title('$\beta_2$','Interpreter','latex')
nexttile(7)
title('$\beta_3$','Interpreter','latex')
nexttile(9)
title('$\beta_4$','Interpreter','latex')
%%
plot_filename='linear_regression_histogram'
ht=2; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

%%

% 
% 
% ss=randsample(1:numel(nreps),1);
% tiledlayout(2,1)
% nexttile
% nonlinfit(zts_unif, Xdat_unif(ss,:), 12, 5);
% nexttile
% nonlinfit(zts_nu, Xdat_nu(ss,:), 12, 5);

res_nu=cell(nreps,1);
res_unif=cell(nreps,1);
tic
parfor ii=1:nreps
    res_nu{ii}=nonlinfit(zts_nu, Xdat_nu(ii,:), 12, 5);
    res_unif{ii}=nonlinfit(zts_unif, Xdat_unif(ii,:), 12, 5);
end
toc

% Extract parameter estimates
for ii=1:numel(res_nu)
    res_nu{ii}=[res_nu{ii}.a0 res_nu{ii}.a1 res_nu{ii}.a2 res_nu{ii}.a3 res_nu{ii}.a4 res_nu{ii}.per1 res_nu{ii}.per2];
end
res_nu=cell2mat(res_nu);
for ii=1:numel(res_unif)
    res_unif{ii}=[res_unif{ii}.a0 res_unif{ii}.a1 res_unif{ii}.a2 res_unif{ii}.a3 res_unif{ii}.a4 res_unif{ii}.per1 res_unif{ii}.per2];
end
res_unif=cell2mat(res_unif);


%% Bottom half plot period estimation part
tiledlayout(2,10,'TileSpacing','tight','Padding','tight')
nexttile(1,[2 3])
min_bin =min(vertcat(res_unif(:,6),res_nu(:,6)));
max_bin =max(vertcat(res_unif(:,6),res_nu(:,6)));
histogram(res_unif(:,6),ceil(sqrt(nreps)),'LineStyle','none',BinLimits=[min_bin max_bin])
hold on
histogram(res_nu(:,6),ceil(sqrt(nreps)),'LineStyle','none',BinLimits=[min_bin max_bin])
xline(avec(6),'LineWidth',1,'color','black')
hold off
ylabel('count','interpreter','latex')
xlabel('NLR estimate','interpreter','latex')
xlim([0 24])
title('$T_1$','Interpreter','latex')

nexttile(4,[2 3])
min_bin =min(vertcat(res_unif(:,6),res_nu(:,7)));
max_bin =max(vertcat(res_unif(:,6),res_nu(:,7)));
set(gca,'YTickLabel',[]);
xlabel('count','interpreter','latex')
histogram(res_unif(:,7),ceil(sqrt(nreps)),'LineStyle','none',BinLimits=[min_bin max_bin])
hold on
histogram(res_nu(:,7),ceil(sqrt(nreps)),'LineStyle','none',BinLimits=[min_bin max_bin])
xlim([0 24])
xline(avec(7),'LineWidth',1,'color','black')
title('$T_2$','Interpreter','latex')
xlabel('NLR estimate','interpreter','latex')

nexttile(7,[1 4])
t=0:.01:24;
ss=randsample(1:numel(nreps),1);
plot(zts_unif,Xdat_unif(ss,:),'.k','color',[0 0.4470 0.7410])
hold on
res=res_unif;
plot(t,res(ss,1)+res(ss,2)*sin(2*pi*t/res(ss,6))+ ...
res(ss,3)*cos(2*pi*t/res(ss,6))+res(ss,4)*sin(2*pi*t/res(ss,7))+ ...
res(ss,5)*cos(2*pi*t/res(ss,7)),'color',[0 0.4470 0.7410])
xlim([0 24])
yticks([-1.5 1.5])
xticks([0 12 24])
ylim([-1.5 1.5])
set(gca,'XTickLabel',[]);
ylabel('$y(t)$','Interpreter','latex')

hold off
nexttile(17,[1 4])
plot(zts_nu,Xdat_nu(ss,:),'.r','color',[0.8500 0.3250 0.0980])
hold on
res=res_nu;
plot(t,res(ss,1)+res(ss,2)*sin(2*pi*t/res(ss,6))+ ...
res(ss,3)*cos(2*pi*t/res(ss,6))+res(ss,4)*sin(2*pi*t/res(ss,7))+ ...
res(ss,5)*cos(2*pi*t/res(ss,7)),'color',[0.8500 0.3250 0.0980])
hold off
xlim([0 24])
xticks([0 12 24])
ylabel('$y(t)$','Interpreter','latex')
yticks([-1.5 1.5])
ylim([-1.5 1.5])
xlabel('$t$','Interpreter','latex')

plot_filename='nonlinear_regression_histogram_top'
ht=2; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

%%
% Second half do histograms
tiledlayout(2,10,'TileSpacing','tight')
for i=1:5
    min_bin =min(vertcat(res_unif(:,i),res_nu(:,i)));
    max_bin =max(vertcat(res_unif(:,i),res_nu(:,i)));
    nexttile(1+2*(i-1),[2 2])
    histogram(res_unif(:,i),ceil(sqrt(nreps)),'LineStyle','none','FaceColor','black','FaceAlpha',.6,BinLimits=[min_bin max_bin])
    hold on
    ylim([0 1000])
    histogram(res_nu(:,i),ceil(sqrt(nreps)),'LineStyle','none','FaceColor',[0 .44 1],'FaceAlpha',.6,BinLimits=[min_bin max_bin])
    xline(avec(i),'LineWidth',1,'color','green')
    if i>1
        set(gca,'YTickLabel',[]);
    end
end


nexttile(1)
ylabel('count')
nexttile(1+2*(3-1))

xlabel('non-linear regression estimate','interpreter','latex')
legend({'uniform','non-uniform','true value'},'Location','southoutside','NumColumns',3)




%%
plot_filename='nonlinear_regression_histogram_bottom'
ht=2; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too


%%
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
