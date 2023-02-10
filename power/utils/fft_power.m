
clear
warning('off','MATLAB:rankDeficientMatrix')
clf

simtype='medium';
checkUseGPU
tic
param.NL=4;
param.NR=4;
Nacro=32;
param.Nacro=Nacro;
param.freq_true=3.9; % freq used in regression model
param.Amp=1;
param.noise2=1/1.5;

% pwr=simulatePWR(param,'uniform');
% acrovec=linspace(0,2*pi,Nacro+1);
% acrovec=acrovec(1:end-1);
% tiledlayout(2,1)
% nexttile
% plot(acrovec,pwr)
% ylim([0,1])
% [m,f,P1]=get_max_psd(pwr)
% nexttile
% plot(f,P1)

noisevals=[2 1.5 1 0.5];
N_SNR=length(noisevals);
Nfreq=2^7;
freqmin=1;
freqmax=16;
freq_vals=linspace(freqmin,freqmax,Nfreq);
cvals=linspace(0.8,0,length(noisevals))'.*[1 1 1];

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


for ii=1:N_SNR
    fpower=NaN(1,Nfreq);
    paramloc=param;
    paramloc.noise2=noisevals(ii);
    for jj=1:Nfreq
        paramloc.freq_true=freq_vals(jj);
        pwr=simulatePWR(paramloc,'uniform');
        fpower(jj)=get_max_psd(pwr);
    end
    plot(freq_vals,fpower,'color',cvals(ii,:))
    hold on
    xlim([1,freqmax])
    ylim([0 0.5])
    drawnow
end

warning('on','MATLAB:rankDeficientMatrix')
ylabel('$\mathrm{max}(|a_{n}|)$','interpreter','latex')
xlabel('$f$','interpreter','latex')



h=colorbar;
ax=gca;
ax.CLim=([min(noisevals),max(noisevals)]);
cvals_plt=linspace(cvals(1,1),cvals(end,1),100)'.*[1 1 1];
colormap(flipud(cvals_plt))
set(h,'TickLabelInterpreter','latex')
h.Label.String='$n$';
h.Label.String='$\sigma$';
h.Label.Interpreter='latex';
h.Ticks=sort(noisevals);

plot_filename=strcat('figs/fft_power_',simtype);

ht=2; % height
wd=5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

function [mpsd,f,P1]=get_max_psd(X)
N=length(X);
Fs = N;             % Sampling frequency                    
T  = 1/Fs;          % Sampling period       
L  = N;             % Length of signal
t  = (0:L-1)*T;     % Time vector

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
mpsd=max(abs(P1(2:end)));
end
