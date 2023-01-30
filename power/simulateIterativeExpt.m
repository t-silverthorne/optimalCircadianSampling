% test if uniform vs non-uniform sampling improves performance when data is
% generated using randomly obtained Fourier coefficients
close all
clear
simtype='medium';
checkUseGPU
param.NL=5;
param.NR=3;
param.freq_true=2; % period used in regression model
param.Amp=1;
param.noise2=1/1.5;

Numdelta=1;
deltavals=0;
for j=1:Numdelta
    nexttile(j)
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
% param.freqtrue=param.freq;
% param.Amp=1.5;