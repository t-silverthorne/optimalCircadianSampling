% path and data
close all
clear all
load('data/sweepNvals_test.mat')
addpath('utils_core')
addpath('utils_cost_fun')

% aesthetics 
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fig_width = 10; % width in inches
fig_height = 5; % height in inches
fig_aspect_ratio = fig_width / fig_height;
fig = figure;
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);


tiledlayout(1,2)
nexttile(1)
[M,c]=contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.5)
set(gca,'XScale','log')



p.Nbatch=1;
p.Nacro=4;
p.Amp=1;
p.freq=3.7;
[pwr,~,~,~]=getPower(t_unif,p)

xline()
nexttile(2)
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1); % get acrophases
plot(acrovec,reshape(pwr,[1 p.Nacro 1 1]),'-k')

drawnow
