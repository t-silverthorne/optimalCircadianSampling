addpath('../utils_core')
addpath('../utils_cost_fun')
clear all
clf
p.Nmeas     = 8;
p.Nresidual = 5e3;
p.Nperm     = 2e2;
p.freq      = 3.7;
p.Amp       = 2;
p.noise     = 1;
acro        = 0; % choose random acrophase
p.Nbatch    = 2;
%[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
tiledlayout(2,1)
nout=50;
nin=10;
for pind=1:nout
t = sort(rand(1,p.Nmeas));

for kk=1:nin
    eps=randn(p.Nresidual,p.Nmeas);
    Y=p.Amp*cos(2*pi*t*p.freq-acro)+p.noise*eps;
    X=constructX(t,p);
    
    betas_obs = (X'*X)\(X'*Y'); %pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); 
    phi_est   = atan2(betas_obs(2,:),betas_obs(3,:));
    amp_est   = sqrt(betas_obs(2,:).^2 + betas_obs(3,:).^2);
    
    fits_obs  = pagetranspose(pagemtimes(X,betas_obs));
    SSres_obs = sqrt(sum((fits_obs-Y).^2,2));
    
    YI = repmat(Y,[1 1 p.Nperm]);
    for jj=1:p.Nresidual
        for ii=1:p.Nperm
            YI(jj,:,ii)=YI(jj,randperm(p.Nmeas,p.Nmeas),ii);
        end
    end
    
    betas   = pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
    fits    = pagetranspose(pagemtimes(X,betas));
    SSres   = sqrt(sum((fits-YI).^2,2));
    bin_mat  = SSres<SSres_obs;
    pwr_est = sum(sum(bin_mat,3)/p.Nperm<.05)/p.Nresidual;

    kk
    amp_stat_vec(kk)  = mean(amp_est-p.Amp);
    phase_stat_vec(kk) =abs(mean(exp(1j*phi_est)));
    pwr_est_vec(kk)   =pwr_est;
end
nexttile(1)
qx=quantile(amp_stat_vec,3);
xv    = mean(amp_stat_vec);
xlower=xv-qx(1);
xupper=qx(3)-xv;


qy=quantile(pwr_est_vec,3);
yv = mean(pwr_est_vec);
ylower=yv-qy(1);
yupper=qy(3)-yv;
errorbar(xv,yv,ylower,yupper,xlower,xupper,'.k')
hold on

nexttile(2)
qx=quantile(phase_stat_vec,3);
xv    = mean(phase_stat_vec);
xlower=xv-qx(1);
xupper=qx(3)-xv;
errorbar(xv,mean(pwr_est_vec),ylower,yupper,xlower,xupper,'.k')
hold on
drawnow
end