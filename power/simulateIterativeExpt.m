% test if uniform vs non-uniform sampling improves performance when data is
% generated using randomly obtained Fourier coefficients
close all
clear
param.NL=5;
param.NR=3;
param.useGPU=true;
param.per=2; % period used in regression model
param.Amp=3;

simtype='fast';
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nresidual=30; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'fast'
        param.Nperm=1e2;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'long'
        param.Nperm=1e3;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'verylong'
        param.Nperm=1e3;
        param.Nresidual=5e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
end

Numdelta=1;
tiledlayout(1,Numdelta)
deltavals=0;%linspace(-5,5,Numdelta);
for j=1:Numdelta
    nexttile(j)
    param.pertrue=(1+deltavals(j)/100)*param.per;

    pwrUnif=simulatePWR(param,'uniform');
    plot(linspace(0,2*pi,param.Nacro),pwrUnif)
    hold on
    
    pwrNU=simulatePWR(param,'non-uniform');
    plot(linspace(0,2*pi,param.Nacro),pwrNU)
    
    legend({'unif','nu'})
    ylim([0,1])
    drawnow
end






% Uniform oscillates, others do not
% param.NL=8; param.NR=4;
% param.per=18; % period used in regression model
% param.pertrue=param.per;
% param.Amp=1.5;