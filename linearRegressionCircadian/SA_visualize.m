close all
f=openfig('SA_good_run_1000.fig')

ax=findobj(f,'type','axes')

h1 = findobj(ax(1),'Type','line')
h2 = findobj(ax(2),'Type','line')

%%
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure(2)
tiledlayout(2,1,'TileSpacing','tight')

nexttile(1)
plot(h2.XData,h2.YData,'.k')    
hold on
ylabel('current $\Psi(\xi)$')
xlim([0 1000])
yline(2.4090,'--k','Color',[0 0.4470 0.7410],'LineWidth',1.3)
hold off
set(gca,'XTickLabel',[]);

nexttile(2)
plot(h1.XData,h1.YData,'.k')    
xlabel('iteration')
ylabel('best $\Psi(\xi)$')
ll=yline(2.4090,'--k','Color',[0 0.4470 0.7410],'LineWidth',1.3,'DisplayName','$M(\xi^N_{\mathrm{unif}}$)')
xlim([0 1000])
legend(ll,'location','northeast','Interpreter','latex')

plot_filename='simulated_annealing_best'
ht=3; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too


