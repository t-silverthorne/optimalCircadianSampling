tic
A1=10;
A2=10;
Nleft=8;
Nright=12;
tau=1/3;
nreps=100;
per1=.4;
per2=1;

Ntimes=Nleft+Nright;
mt1=linspace(0,tau,Nleft+1);
mt1=mt1(1:end-1);
mt2=linspace(tau,1,Nright+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2]; % construct non-uniform grid 
mt_unif=linspace(0,1,Ntimes+1); % construct uniform grid 
mt_unif=mt_unif(1:end-1);

dphi1=2*pi/(10+1);
dphi2=dphi1;
phi1vals=0:dphi1:2*pi;
phi2vals=0:dphi2:2*pi;

phi1=0;
phi2=0.6;

linParams=get_linear_params([A1 A2 phi1 phi2]);

a1=linParams(:,1);
a2=linParams(:,2);
a3=linParams(:,3);
a4=linParams(:,4);
get_Xdat = @(zts) a1.*sin(2*pi*zts/per1)+a2.*cos(2*pi*zts/per1) + ...
    a3.*sin(2*pi*zts/per2)+a4.*cos(2*pi*zts/per2) + randn(nreps,numel(zts));

nu_noisy=true;
unif_noisy=true;
period_er_nu=[];
period_er_unif=[];
niter=0;
while (nu_noisy || unif_noisy) && niter < 10
    niter=niter+1;
    Xdat_unif=get_Xdat(mt_unif);
    Xdat_nu=get_Xdat(mt_nu);
    
    if nu_noisy
        [period_er_nu_loc]  =get_mean_period_error(mt_nu,Xdat_nu,per1,per2);
    end
    if unif_noisy
        [period_er_unif_loc]=get_mean_period_error(mt_unif,Xdat_unif,per1,per2);
    end

    if numel(period_er_nu)>0
        period_er_nu=[period_er_nu period_er_nu_loc];
        period_er_unif=[period_er_unif period_er_unif_loc];
    else
        period_er_nu=period_er_nu_loc;
        period_er_unif=period_er_unif_loc;
    end
    
    mean_er_nu=mean(period_er_nu);
    mean_er_unif=mean(period_er_unif);
    snorm_nu=std(period_er_nu)/mean(mean_er_nu);
    snorm_unif=std(period_er_unif)/mean(mean_er_unif);
    fprintf('unif noise:  %d\n',snorm_unif);
    fprintf('unif    er:  %d\n',mean_er_unif);
    fprintf('nu   noise:  %d\n',snorm_nu);
    fprintf('nu      er:  %d\n\n',mean_er_nu);
    unif_noisy= (snorm_unif>1e-1);
    nu_noisy= (snorm_nu>1e-1); 
end
toc