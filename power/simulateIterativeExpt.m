% test if uniform vs non-uniform sampling improves performance when data is
% generated using randomly obtained Fourier coefficients
close all
clear
param.NL=5;
param.NR=3;
param.useGPU=false;
param.freq=1; % period used in regression model
param.Amp=1;

simtype='rough';
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nresidual=30; 
        param.Nacro=32;
    case 'fast'
        param.Nperm=1e2;
        param.Nresidual=1e2;
        param.Nacro=32; 
    case 'long'
        param.Nperm=1e3;
        param.Nresidual=1e3; 
        param.Nacro=32; 
    case 'verylong'
        param.Nperm=1e3;
        param.Nresidual=5e3; 
        param.Nacro=32; 
end

max(abs(simulatePWR(param,'uniform')))


% Numdelta=1;
% deltavals=0;
% for j=1:Numdelta
%     nexttile(j)
%     param.per=(1+deltavals(j)/100)*param.per;
% 
%     pwrUnif=simulatePWR(param,'uniform');
%     plot(linspace(0,2*pi,param.Nacro),pwrUnif)
%     hold on
%     
%     pwrNU=simulatePWR(param,'non-uniform');
%     plot(linspace(0,2*pi,param.Nacro),pwrNU)
%     
%     legend({'unif','nu'})
%     ylim([0,1])
%     drawnow
% end






% Uniform oscillates, others do not
% param.NL=8; param.NR=4;
% param.per=18; % period used in regression model
% param.freqtrue=param.freq;
% param.Amp=1.5;