%% preamble
close all
c1="#2271B2"; % blue uniform colour
c2="#F748A5"; % pink 2-uniform colour
c3="#D55E00"; % orange Unif(0,1) colour
c4="#359B73"; % green jittered colour
addpath('../utils')
fig=figure;
tiledlayout(4,4,'TileSpacing','compact')
%%
p.Amp       = 2;
p.freq      = 3.6;
p.Nmeas     = 8; 
p.mesor     =0
p.Nresidual = 1e3;
    
%nout           = [50 50 50 1 1];
nout           = [50 50 50 50 1];
designs_to_run = {'random','2-uniform','jittered','optimal','uniform'};

% nout           = [10 10 1];
% designs_to_run = {'2-uniform','optimal','uniform'};




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
                options=optimoptions('fmincon',...
                  'CheckGradients',false,...
                  'SpecifyObjectiveGradient',true);
                [t,~]=fmincon(@(t) get_MINUS_MinLambdaAndDiff(t,p.freq),rand(1,p.Nmeas),[],[],[],[],zeros(1,p.Nmeas),ones(1,p.Nmeas),[],options);
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
        Jvec = get_Jvec(t,p);
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

%clean up aesthetics
linkaxes([nexttile(5) nexttile(6)],'y')
linkaxes([nexttile(9) nexttile(10) nexttile(11)],'y')
linkaxes([nexttile(13) nexttile(14) nexttile(15) nexttile(16)],'y')

linkaxes([nexttile(1) nexttile(5) nexttile(9) nexttile(13)],'x')
linkaxes([nexttile(6) nexttile(10) nexttile(14)],'x')
linkaxes([nexttile(11) nexttile(15)],'x')

leg=legend('Location','southoutside','NumColumns',3)
leg.Layout.Tile='south'
set(findall(gcf,'-property','FontSize'),'FontSize',10)

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

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT for paper
plot_filename='/Users/turnersilverthorn/research/overleaf/samplingPaper/figures/fig3.png';
ht=4; % height
wd=5.8; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename, ...
    '-dpng','-r600') % -r sets the resolution

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT for presentation
fig_width  = 4; % width in inches
fig_height = 2.9; % height in inches
plot_filename='/Users/turnersilverthorn/research/overleaf/samplingTalk/figsAMMCS/fig3.png';
set(fig, 'PaperUnits', 'inches', ...
         'PaperSize', [fig_width fig_height], ...
         'PaperPosition', [0 0 fig_width fig_height], ...
         'PaperPositionMode', 'manual', ...
         'Units', 'inches', ...
         'Position', [0 0 fig_width fig_height]);
set(findall(gcf,'-property','FontSize'),'FontSize',6)
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution



function Jvec = get_Jvec(mt,p)
    N=length(mt);
    Amp  = p.Amp;
    freq = p.freq;
    beta1= p.mesor;
    numphi=32;
    phis = linspace(0,2*pi,numphi);
    for ii=1:numphi
        phi  = phis(ii);
        beta2=Amp*cos(phi);
        beta3=Amp*sin(phi);
        fit=fit_cosinor_model(beta1+ beta2*cos(2*pi*freq*mt)+beta3*sin(2*pi*freq*mt)+randn(p.Nresidual,N),mt,1/freq);
        Jvec(ii,1)  = mean(fit.amplitudes-Amp); % amp bias
        Jvec(ii,2)  = mean(abs(exp(1j*fit.acrophases_rad) -exp(1j*phi))); % acro bias
        Jvec(ii,3)  = var(fit.amplitudes); % amplitude variance
        Jvec(ii,4)  = 1-abs(mean(exp(1j*fit.acrophases_rad))); % acro variance
        Jvec(ii,5)  = 1-getMinPower(Amp,freq,mt); % worst case power
    end
    Jvec=max(Jvec,[],1);
end