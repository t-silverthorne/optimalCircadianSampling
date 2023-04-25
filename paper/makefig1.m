plot_filename='fig1';
addpath('utils_core')
addpath('utils_cost_fun')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

font_size = 10;
set(0, 'DefaultAxesFontSize', font_size);
set(0, 'DefaultTextFontSize', font_size);

% show grids
clf
ny1=4;ny2=3;
nx1=3;nx2=2;nx3=2;
tileind=@(row,col) sub2ind([nx1+nx2+nx3,ny1+ny2],col,row);

tiledlayout(ny1+ny2,nx1+nx2+nx3,'TileSpacing','tight','padding','none')

%% Show grids
nexttile(tileind(1,1),[ny1,nx1])
addpath('utils_core')
rng(2345)
p.Nmeas=8;
N=p.Nmeas;
p.permMethod='fast';
tauL=0.4;
tauR=0.6;
[t_unif,t_nu]=getSamplingSchedules(4,4,tauL,tauR);
t_rand=sort(rand(1,N));

t_jit = t_unif + 1e-2*rand(1,N);
c1="#648FFF"; % uniform colour
c2="#785EF0"; % 2-uniform colour
c3="#DC267F"; % Unif(0,1) colour
c4="#FE6100"; % jittered colour
plot(t_unif,4,'ok','MarkerEdgeColor','none','MarkerFaceColor',c1)

hold on
plot(t_nu,3,'ok','MarkerEdgeColor','none','MarkerFaceColor',c2)
plot(t_rand,2,'ok','MarkerEdgeColor','none','MarkerFaceColor',c3)
plot(t_jit,1,'ok','MarkerEdgeColor','none','MarkerFaceColor',c4)

yv=.2*linspace(-1,1,10);

plot(repmat(tauL,1,10),3+yv,'-k','Color',c2)
plot(repmat(tauR,1,10),3+yv,'-k','Color',c2)
ylim([0.5 4.5])
xlabel('$t$','interpreter','latex')
yticks([1 2 3 4])
yticklabels({'jittered','random','2-uniform','uniform'})

p.Nacro     = 64;
p.freq      = 3.8;
p.Amp       = 1.5;
p.noise     = .5;
p.Nbatch    = 20;
p.Nresidual = 1;
p.Nperm     = 1e2;

t_hires=linspace(0,1,100);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);


Xhr=constructX(t_hires,p);
acind=4;

%% Show simulated data and fits
for ii=1:4
    nexttile(tileind(ii,nx1+1),[1,nx2])
    switch ii
        case 1
            tloc=t_unif;
            cloc=c1;
        case 2
            tloc=t_nu;
            cloc=c2;
        case 3
            tloc=t_rand;
            cloc=c3;
        case 4
            tloc=t_jit;
            cloc=c4;
    end
    [beta,Y]=get_beta(tloc,p);
    plot(tloc,Y(:,:,:,acind),'ok','MarkerEdgeColor','none','MarkerFaceColor',cloc)
    hold on
    plot(t_hires,Xhr*beta(:,:,:,acind),'-k','color',cloc)
    ylim([-3 3])
    if ii==4
        xlabel('$t$','Interpreter','latex')
    else
        set(gca,'XTickLabel',[]);
    end
    ylabel('$y(t)$','Interpreter','latex')
end
p.Nresidual = 1e3;
%% Plot power and counts
[I3,I4]=constructUtilMats(p);
for ii=1:4
    nexttile(tileind(ii,nx1+nx2+1),[1,nx3])
    switch ii
        case 1
            tloc=t_unif;
            cloc=c1;
        case 2
            tloc=t_nu;
            cloc=c2;
        case 3
            tloc=t_rand;
            cloc=c3;
        case 4
            tloc=t_jit;
            cloc=c4;
    end

    [pwr_est,amp_est,phi_est] = getPowerBatch(tloc,p,I3,I4);
    histogram(amp_est(:,acind),'EdgeColor','none','FaceColor',cloc,'Normalization','probability')

    hold on
    xline(p.Amp,'--k')
    xlim([0,5])
    if ii==4
        xlabel('$\hat{A}$')
    else
        set(gca,'XTickLabel',[]);
    end
    ylabel('$p(\hat{A})$','interpreter','latex')

    

    nexttile(tileind(ny1+1,1),[ny2,nx1])
    acrovec=linspace(0,2*pi,p.Nacro+1);
    acrovec=acrovec(1:end-1); % get acrophases
    plot(acrovec,reshape(mean(pwr_est,1),1,p.Nacro),'-k','Color',cloc)
    ylim([0,1])
    xlim([0,2*pi])
    % x axis for angular variable
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})
    xlabel('$\phi$','interpreter','latex')
    ylabel('power $\gamma(\phi)$','interpreter','latex')

    hold on
    drawnow
end
%%
nexttile(tileind(ny1+1,nx1+1),[ny2,nx2+nx3])
p.Nresidual = 1;
p.Nperm     = 1;
p.Nbatch    = 1;
[I3,I4]=constructUtilMats(p);
%%
for ii=1:4
    ploc=p;
    switch ii
        case 1
            tloc=t_unif;
            cloc=c1;
        case 2
            tloc=t_nu;
            cloc=c2;
        case 3
            tloc=t_rand;
            cloc=c3;
        case 4
            tloc=t_jit;
            cloc=c4;
        case 5 
            tloc=0:.001:1;
            ploc.Nmeas=length(tloc);
            cloc='black';
    end
    
    Y=getSimulatedData(tloc,ploc);
    Y=Y(:,:,:,acind);
    %[pxx,f]=periodogram(Y',[],[],ploc.Nmeas);
    [~,ind]=sort(tloc);
    [pxx,f]=plomb(Y(ind)',tloc(ind),[],20);
    if ii<5
        plot(f,pxx,'-k','HandleVisibility','off','color',cloc)
    else
        plot(f,pxx,'-k','HandleVisibility','off','color',cloc,'LineWidth',2)
    end
    hold on
    if ii==4
        xline(p.freq,'--k')
    end
    drawnow
    xlim([0,4])
    xlabel('$f$','interpreter','latex')
    ylabel('spectral density $\rho(f)$','interpreter','latex')
end
%fits_obs  = pagetranspose(pagemtimes(X,betas_obs));
%SSres_obs = sqrt(sum((fits_obs-Y).^2,2));
%%
plot_filename='figs/fig1'
ht=8; % height
wd=12; % width
% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
fig=gcf;ax=fig.CurrentAxes;fig.Color='w';fig.OuterPosition=fig.InnerPosition;
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

function [beta,Y]=get_beta(t,p)
Y=getSimulatedData(t,p);
X=constructX(t,p);
beta=pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); % observed error
end