plot_filename='fig1';
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% show grids
clf
ny1=4;ny2=4;
nx1=3;nx2=2;nx3=1;
tileind=@(row,col) sub2ind([nx1+nx2+nx3,ny1+ny2],col,row);

tiledlayout(ny1+ny2,nx1+nx2+nx3,'TileSpacing','compact')
nexttile(tileind(1,1),[ny1,nx1])
addpath('utils_core')
rng(2345)
p.Nmeas=8;
N=p.Nmeas;
tauL=0.4;
tauR=0.6;
[t_unif,t_nu]=getSamplingSchedules(4,4,tauL,tauR);
t_rand=sort(rand(1,N));

t_jit = t_unif + 5e-2*rand(1,N);
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
xlabel('time','interpreter','latex')
yticks([1 2 3 4])
yticklabels({'jittered','uniformly random','2-uniform','uniform'})
% ax1=gca;
% figure
% tiledlayout(2,1)
% t1=nexttile(1)
% copyobj(ax1,t1)
% t2=nexttile(2)
% copyobj(get(ax1,'children'),t2)
% difference in sampling


% p.Nacro     = 16;
% p.Nresidual = 1;
% p.freq      = 5.1;
% p.Amp       = 1;
% p.noise     = 0.5;
% p.Nbatch    = 1;

p.Nacro     = 16;
p.Nresidual = 1;
p.freq      = 3.8;
p.Amp       = 1.5;
p.noise     = .5;
p.Nbatch    = 1;

t_hires=linspace(0,1,100);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);


Xhr=constructX(t_hires,p);
acind=4;
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
end
p.Nresidual = 5e3;
p.Nperm     = 1e2;
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
    [pwr_est,amp_est,phi_est] = getPower(tloc,p,I3,I4);
    histogram(amp_est(:,:,:,acind),'EdgeColor','none','FaceColor',cloc,'Normalization','probability')
    hold on
    xline(p.Amp,'--k')
    xlim([0,5])
    

    nexttile(tileind(ny1+1,1),[ny1,nx1])
    acrovec=linspace(0,2*pi,p.Nacro+1);
    acrovec=acrovec(1:end-1); % get acrophases
    plot(acrovec,reshape(pwr_est,1,p.Nacro),'-k','Color',cloc)
    ylim([0,1])
    hold on
    drawnow
end
%%
nexttile(tileind(ny1+1,nx1+1),[ny2,nx2+nx3])
p.Nresidual = 1;
p.Nperm     = 1;
[I3,I4]=constructUtilMats(p);

for ii=1:5
    switch ii
        case 1
            tloc=t_unif;
            ploc=p;
            cloc=c1;
        case 2
            tloc=t_nu;
            ploc=p;
            cloc=c2;
        case 3
            tloc=t_rand;
            ploc=p;
            cloc=c3;
        case 4
            tloc=t_jit;
            ploc=p;
            cloc=c4;
        case 5 
            tloc=t_hires;
            ploc=p;
            ploc.Nmeas=length(tloc);
            cloc='black';
    end
    
    Y=getSimulatedData(tloc,ploc);
    Y=Y(:,:,:,acind);
    [pxx,f]=periodogram(Y',[],[],ploc.Nmeas);
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
end
%fits_obs  = pagetranspose(pagemtimes(X,betas_obs));
%SSres_obs = sqrt(sum((fits_obs-Y).^2,2));

function [beta,Y]=get_beta(t,p)
Y=getSimulatedData(t,p);
X=constructX(t,p);
beta=pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); % observed error
end