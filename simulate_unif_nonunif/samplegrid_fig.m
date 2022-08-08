
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% non-uniform grid
Ntimes=31; % Ideal number of samples
frq1=2; % samples per hour in first region
frq2=.5;
f= [frq1 frq2 -1];
A=-[frq1 frq2  0; zeros(2,3)];
b=-[Ntimes-2; zeros(2,1)];
Aeq=[1 1 0; 0 0 1; 0 0 0];
beq=[24; Ntimes;0];
lb=[1 1 1];
nui=intlinprog(f,1:3,A,b,Aeq,beq,lb);
Ntimes=ceil(nui(1:2)'*[frq1 frq2 ]'); % number of samples consistent with restrictions

% finish non-uniform grid
mt1=linspace(0,nui(1),frq1*nui(1)+1);
mt1=mt1(1:end-1);
mt2=linspace(nui(1),nui(1)+nui(2),ceil(frq2*nui(2))+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2];

% uniform grid 
mt_unif=linspace(0,24,Ntimes+1);
mt_unif=mt_unif(1:end-1);


clf
tiledlayout(1,1,'TileSpacing','none')
nexttile
scatter(mt_unif,1,'.k')
hold on
scatter(mt_nu,0,'.k')
ylim([-0.5 1.5])
xlim([0 24])
xticks([0 6 12 18 24])
xline(12,'--k')
xlabel('hours')
yticks([0 1])
yticklabels({'uniform','non-uniform'})

plot_filename='samplegrid_fig'
ht=1.5; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too











