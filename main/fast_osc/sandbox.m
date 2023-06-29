addpath('../utils/')
Nmeas=8;
%% power as a function of acrophase
[mt,~]   = getSamplingSchedules(Nmeas,0,0,0);

acros   = linspace(0,2*pi,2^7);
powvals = arrayfun(@(x) getPowerExact(2,x,3.8,mt),acros);
plot(acros,powvals,'-ok')


%% 2D grid
tic
numFgrid=2^7;
numAgrid=2^4;
freqvals        = linspace(1,18,numFgrid);
Ampvals         = logspace(-1,1,numAgrid);

[Agr,Fgr]=meshgrid(Ampvals,freqvals);
Pmin=arrayfun(@(Agr,Fgr) getMinPower(Agr,Fgr,mt),Agr, Fgr);
[M,c]=contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.5)
set(gca,'XScale','log')
xlabel('amplitude')
toc


function beta=getPowerExact(Amp,acro,freq,mt)
Nmeas    = length(mt); % num samples
alpha    = 0.05; % type I error rate

csq    = cos(2*pi*freq*mt-acro);
lambda   =  Amp^2*csq*csq'; % non-centrality

beta=1-ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);
end

function betamin=getMinPower(Amp,freq,mt)
[~,betamin]=fminbnd(@(phi) getPowerExact(Amp,phi,freq,mt),0,2*pi);
end