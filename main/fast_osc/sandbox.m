close all
clear
addpath('../utils/')
% settings
Nmeas=8;
Amp=2;
freq=3;
maxit=1e3;

fast_run=true;

% power as a function of acrophase
[mt,~]   = getSamplingSchedules(Nmeas,0,0,0);

acros   = linspace(0,2*pi,2^5);
powvals = arrayfun(@(x) getPowerExact(Amp,x,freq,mt),acros);
plot(acros,powvals,'-ok')
hold on

% optimize power
opts = optimoptions(@simulannealbnd,'Display','iter', ...
            'MaxIterations',maxit,'DisplayInterval',100,'ReannealInterval',50);
x0=rand(1,Nmeas);
xopt = simulannealbnd(@(t) -getMinPower(Amp,freq,t), ...
		       x0,zeros(1,Nmeas),ones(1,Nmeas),opts);

@(t) getMinPower(Amp,freq,t)

acros   = linspace(0,2*pi,2^5);
powvals = arrayfun(@(x) getPowerExact(2,x,3.8,xopt),acros);
plot(acros,powvals,'-ob')
xlim([0,2*pi])
ylim([0,1])
%% 2D grid
tic
numFgrid=2^7;
numAgrid=2^4;
if fast_run
    numFgrid=2^5;
    numAgrid=2^3;
end
freqvals        = linspace(1,18,numFgrid);
Ampvals         = logspace(-1,1,numAgrid);

[Agr,Fgr]=meshgrid(Ampvals,freqvals);
Pmin=arrayfun(@(Agr,Fgr) getMinPower(Agr,Fgr,mt),Agr, Fgr);
[M,c]=contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.5)
set(gca,'XScale','log')
xlabel('amplitude')
toc




