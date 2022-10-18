plot_filename='bayesian1Dexample';
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%%
tiledlayout(1,3,'TileSpacing','compact')
first=true;
for jj=[1/4 1/3 1/2]
    nexttile
    for i=3:10
        [unif_vals,sd_unif,nu_vals,sd_nu]=estimateBayesianFIMdetCirc(i,2+i,0,jj,'sdev');
        mean_unif=mean(unif_vals);
        mean_nu=mean(nu_vals);
        errorbar(i,mean_unif,sd_unif,'.k','MarkerSize',5)
        hold on
        errorbar(i,mean_nu,sd_nu,'.b','MarkerSize',5)
        drawnow
    end
    xlabel('$N_L$','interpreter','latex')
    switch jj
        case 1/4
            ylabel('$S(\xi_{[0,\tau]}(N_L,N_L+2))$','interpreter','latex')
            title('($\tau=\frac{1}{4}$)','interpreter','latex')
            xticks(3:2:10)
        case 1/3
            title('($\tau=\frac{1}{3}$)','interpreter','latex')
            set(gca,'YTickLabel',[]);
            legend({'uniform','non-uniform'},'Location','southoutside','NumColumns',2)
            xticks(3:2:10)
        case 1/2
            title('($\tau=\frac{1}{2})$','interpreter','latex')
            set(gca,'YTickLabel',[]);
            xticks(3:2:10)
            
    end
    hold off
    ylim([8 20])
end

ht=2.5; % height
wd=6.5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

%% heatmap part
close all
tiledlayout(1,3,'TileSpacing','tight')
nii=6;
njj=6;
unif=NaN(nii,njj);
nu=NaN(nii,njj);
for taulim=[1/4 1/3 1/2]
    for ii=1:nii
        disp(ii)
        for jj=1:njj
            [unif_vals,~,nu_vals,~]=estimateBayesianFIMdetCirc(3+ii,3+jj,0,taulim,'sdev');
            unif(ii,jj)=mean(unif_vals);
            nu(ii,jj)=mean(nu_vals);
        end
    end
    nexttile
    h=heatmap(3+(1:njj),3+(nii:-1:1),100*flipud((real(nu)-real(unif))./real(unif)));
    h.GridVisible = 'off';
    switch taulim
        case 1/4
            h.ColorbarVisible = 'off';
            h.XLabel='$N_R$';
            h.YLabel='$N_L$';            
        case 1/3
            h.ColorbarVisible = 'off';
            h.XLabel='$N_R$';
            Ax = gca;
            Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
        case 1/2
            h.XLabel='$N_R$';
            h.NodeChildren(2).TickLabelInterpreter='latex';
            h.NodeChildren(2).Label.String='$S(\xi_{[0,\tau]}(N_L,N_R))$';
            h.NodeChildren(2).Label.Interpreter='latex';
            Ax = gca;
            Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
            h.ColorbarVisible
    end
    h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
    h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';

    drawnow
    h.ColorLimits = [-100 100];
end
%%
nexttile(2)
h=heatmap(1:njj,1:nii,real(nu))
h.YDisplayLabels=''
h.GridVisible = 'off';
h.ColorLimits = [10 50];
