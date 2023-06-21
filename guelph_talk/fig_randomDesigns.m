%% preamble
load('data/guelph_opt.mat')
close all
testing=false;
def_colours
addpath('utils_core')
addpath('utils_cost_fun')

%% aesthetics 
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fig_width  = 5; % width in inches
fig_height = 3.25; % height in inches
fig_aspect_ratio = fig_width / fig_height;
fig = figure;
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);
tiledlayout(4,4,'TileSpacing','compact')
set(findall(gcf,'-property','FontSize'),'FontSize',4)

%% main
if testing
    p.Nbatch    = 1;
    p.Nperm     = 1e1;
    p.Nresidual = 1e1;
    nout           = [3 3 3 3 1];
    designs_to_run = {'random','2-uniform','jittered','optimal','uniform'};
else
    p.Nbatch    = 1;
    p.Nperm     = 1e2;
    p.Nresidual = 1e3;
    p.permMethod       = 'FY_double_for_loop'; %p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
    p.permActionMethod = 'index';
    
%     fname='test_rand';    
%     nout           = [50 50 50 1];
%     designs_to_run = {'random','2-uniform','jittered','uniform'};

    fname='optim_result';
    nout           = [50 50 50 size(xmaster,1) 1];
    designs_to_run = {'random','2-uniform','jittered','optimal','uniform'};
end



lab_list = {'amp bias','acro bias','amp var','acro var', '$1-\textrm{power}$'};

for kk=1:length(designs_to_run)
    for pind=1:nout(kk)
        marktype  = '.';
        switch designs_to_run{kk} % change color/marker size for each design
            case 'uniform'
                [t,~]     = getSamplingSchedules(p.Nmeas,0,0,0);
                cloc      = c1;
                marksize  = 6;
                dname     = 'uniform';
                marktype  = 'pentagram';
            case 'optimal'
                t         = xmaster(pind,:);
                cloc      = 'black';
                marksize  = 6;
                dname     = 'optimal';
            case '2-uniform'
                [~,t]     = getSamplingSchedules(6,2,0,0.3);
                cloc      = c2;
                marksize  = 6;
                dname     = '2-uniform';
            case 'random'
                t         = sort(rand(1,p.Nmeas));
                cloc      = c3;
                marksize  = 6;
                dname     = 'random';
            case 'jittered'
                [t_unif,~]= getSamplingSchedules(p.Nmeas,0,0,0);
                t         = t_unif + 2.5e-2*rand(1,p.Nmeas);
                cloc      = c4;
                marksize  = 6;
                dname     = 'jittered';                
        end
        hvbool='off'; % decide if it should be in legend
        if pind==1
            hvbool='on';
        end
        Jvec = wrap_getCostFun(t,p,active_inds);
        for ii=1:4
            for jj=1:ii
                nexttile(jj + 4*(ii-1))
                if strcmp(designs_to_run{kk},'uniform')
                    cloc_out=[177, 3, 252]./255;
                else
                    cloc_out=cloc;
                end
                loglog(Jvec(jj),Jvec(ii+1),'.k','Color',cloc_out,'MarkerFaceColor',cloc, ...
                    'MarkerSize',marksize, ...
                    'DisplayName',dname,'HandleVisibility',hvbool, ...
                    'Marker',marktype);
                hold on
                xlabel(lab_list{jj})
                ylabel(lab_list{ii+1})
                drawnow
            end
        end
    end
end


linkaxes([nexttile(5) nexttile(6)],'y')
linkaxes([nexttile(9) nexttile(10) nexttile(11)],'y')
linkaxes([nexttile(13) nexttile(14) nexttile(15) nexttile(16)],'y')

linkaxes([nexttile(1) nexttile(5) nexttile(9) nexttile(13)],'x')
linkaxes([nexttile(6) nexttile(10) nexttile(14)],'x')
linkaxes([nexttile(11) nexttile(15)],'x')
%%
leg=legend('Location','southoutside','NumColumns',length(designs_to_run))
leg.Layout.Tile='south'
set(findall(gcf,'-property','FontSize'),'FontSize',8)

for ii=1:4
    for jj=1:ii
        nexttile(jj + 4*(ii-1))
        xl=xlim;
        xticks([xl(1) xl(end)])
        tix=get(gca,'xtick')';
        set(gca,'xticklabel',num2str(tix,'%0.1f'))
    end
end

for ii=1:4
    for jj=1:ii
        nexttile(jj + 4*(ii-1))
        yl=ylim;
        yticks([yl(1) yl(end)])
        tix=get(gca,'ytick')';
        set(gca,'yticklabel',num2str(tix,'%0.1f'))
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',6)
        a = get(gca,'YTickLabel');  
        set(gca,'YTickLabel',a,'fontsize',6)
    end
end

xklist=[1 5 9 6 10 11];
for ii=1:length(xklist)
    nexttile(xklist(ii))
    xlabel('')
    set(gca,'XTickLabel',[])
end

yklist=[6 10 14 11 15 16];
for ii=1:length(yklist)
    nexttile(yklist(ii))
    ylabel('')
    set(gca,'YTickLabel',[])
end





%%
print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/',fname),'-dpng','-r600') % -r sets the resolution
savefig(fname)