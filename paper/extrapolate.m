p.Nmeas     = 8;
p.Nacro     = 1;
p.Nresidual = 1e3;
p.Nperm     = 99;
p.noise     = 1;
p.Nbatch    = 2;
p.Amp       = 2;
p.freq      = 2.8;
[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[I3,I4]=constructUtilMats(p);

Y=getSimulatedData(t,p);
addpath('utils_core/')

% Slow method
p.Nperm=1e4-1;
X=constructX(t,p);
tic
[pwr,~]=getPower(Y,X,p)
toc
% Extrapolating method
%%
Innervec=[2e2-1,100-1,59,39,19];
p.Nperm=Innervec(1);
p.Nresidual=1e3;
Y=getSimulatedData(t,p);
tic
[pwrall,binmat]=getPower(Y,X,p);
binvec=sum(binmat,3);
pvec=NaN(length(Innervec),1);
pvec(1)=pwrall;
for ii=2:length(Innervec)
    pvec(ii)=mean(arrayfun(@(ind) hygecdf(Innervec(ii)*.05,Innervec(1),binvec(ind),Innervec(ii)),1:p.Nresidual));
end

Xq= [ones(1,length(Innervec)); 1./Innervec;  1./Innervec.^2;]';

beta=(Xq'*Xq)\(Xq'*pvec);

beta(1)
toc
function [pwr_est,binmat]=getPower(Y,X,p)
betas_obs = (X'*X)\(X'*Y'); % observed error

fits_obs  = (X*betas_obs)';
SSres_obs = sqrt(sum((fits_obs-Y).^2,2));

YI=NaN(p.Nresidual,p.Nmeas,p.Nperm);

for ii=1:p.Nresidual
    Yloc=Y(ii,:);
    for jj=1:p.Nperm
        YI(ii,:,jj)=Yloc(randperm(p.Nmeas,p.Nmeas));
    end
end

betas   = pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits    = pagetranspose(pagemtimes(X,betas));
SSres   = sqrt(sum((fits-YI).^2,2));
binmat  = SSres<SSres_obs;
pwr_est = sum(sum(binmat,3)/p.Nperm<.05)/p.Nresidual;
end