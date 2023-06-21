% compare results of nonlinear regression with multistart

clear
clf

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

tiledlayout(2,4,'TileSpacing','tight','padding','compact')
plot_option='permin';
other_options='';
mkplots(plot_option,other_options)
plot_option='sse';
mkplots(plot_option,other_options)

lg=legend({'uniform','non-uniform'},'NumColumns',2)
lg.Layout.Tile='South'
%%

%%
set(findall(gcf,'-property','FontSize'),'FontSize',10)

plot_filename='figures/hist_sim_result'
% save
ht=3.8; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too



function mkplots(plot_option,other_options)
fname_list={...
        'final_method2_regularize_1_Nsamp_1_nreps_1000_per1_12_per2_4_Nleft_10_Nright_4.mat', ...
        'final_method2_regularize_1_Nsamp_1_nreps_1000_per1_12_per2_4_Nleft_12_Nright_6.mat', ...
        'final_method2_regularize_1_Nsamp_1_nreps_1000_per1_12_per2_4_Nleft_14_Nright_8.mat', ...
        'final_method2_regularize_1_Nsamp_1_nreps_1000_per1_12_per2_4_Nleft_16_Nright_10.mat'};
for ii=1:4
    nexttile
    X=open(fname_list{ii});
    
    sse_unif=NaN(numel(X.gof_unif),1);
    sse_nu=NaN(numel(X.gof_nu),1);
    
    for i=1:numel(sse_unif)
        sse_unif(i)=X.gof_unif{i}.sse;
        sse_nu(i)=X.gof_nu{i}.sse;
    end
    

    per_unif_min=NaN(numel(X.gof_unif),1);
    per_nu_min=NaN(numel(X.gof_nu),1);
    
    per_unif_max=NaN(numel(X.gof_unif),1);
    per_nu_max=NaN(numel(X.gof_nu),1);
    

    for i=1:numel(per_unif_min)
        per_unif_min(i)=min(X.res_unif{i}.per1,X.res_unif{i}.per2);
        per_nu_min(i)=min(X.res_nu{i}.per1,X.res_nu{i}.per2);

        per_unif_max(i)=max(X.res_unif{i}.per1,X.res_unif{i}.per2);
        per_nu_max(i)=max(X.res_nu{i}.per1,X.res_nu{i}.per2);
    end
    
    switch plot_option
        case 'sse'
            histogram(sse_unif,floor(sqrt(X.nreps)),'BinLimits', ...
                [0 0.5])
            hold on
            histogram(sse_nu,floor(sqrt(X.nreps)),'BinLimits', ...
                [0 0.5])
            hold off
            xlabel('SSE')            
            if ii==1
                ylabel('count')
            end
            ylim([0 200])
            xlim([0 0.5])
            xticks([0 0.2 0.4])
            if ii>1
                set(gca,'YTickLabel',[]);
            end
        case 'permin'
            histogram(per_unif_min,floor(sqrt(X.nreps)),'BinLimits', ...
                [0 4.3])
            hold on
            histogram(per_nu_min,floor(sqrt(X.nreps)),'BinLimits', ...
                [0 4.3])
            if ii==1
                
                ylabel('count')
            end
            hold off
            ylim([0 1000])
            xlim([0 4.3])
            xticks([0 2 4])
            if ii>1
                set(gca,'YTickLabel',[]);
            end
            xlabel('$\textrm{min}(\hat{T}_1,\hat{T}_2)$')
            if ii==1
                title('($N_L=10$, $N_R=4$)')
            elseif ii==2
                title('($N_L=12$, $N_R=6$)')
            elseif ii==3
                title('($N_L=14$, $N_R=8$)')
            elseif ii==4
                title('($N_L=16$, $N_R=10$)')
            end
            
        case 'permax'
            histogram(per_unif_max,floor(sqrt(X.nreps)),'BinLimits', ...
                [min(vertcat(per_nu_max,per_unif_max))   max(vertcat(per_nu_max,per_unif_max))])
            hold on
            histogram(per_nu_max,floor(sqrt(X.nreps)),'BinLimits', ...
                [min(vertcat(per_nu_max,per_unif_max))   max(vertcat(per_nu_max,per_unif_max))])
            if ii==1
                legend({'unif','nu'})
            end
            hold off
            
    end

end

end
