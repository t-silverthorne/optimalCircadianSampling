close all
clear

% load('data/test_opt_full.mat')
% fname='test_opt_full';

load('data/test_opt_full.mat')
fname='test_rand';
%%
testing=false;

addpath('utils_core')
addpath('utils_cost_fun')

% aesthetics 
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


lab_list = {'amp bias','acro bias','amp var','1-acro var', '1-power'};

tiledlayout(4,4,'TileSpacing','compact')
set(findall(gcf,'-property','FontSize'),'FontSize',4)

% p.Nmeas     = 8;
% f1=3.8;
% f2=4.8;
% p.Amp       = 3.5;
nout(1)=1;
nout(2)=size(xmaster,1);
nout(3)=20;
nin = 1;
c1="#648FFF"; % uniform colour
c2="#785EF0"; % 2-uniform colour
c3="#DC267F"; % Unif(0,1) colour
c4="#FE6100"; % jittered colour
c5="#000000";   
cols=[c1,c5,c3];

if testing
    p.Nbatch    = 1;
    p.Nperm     = 1e1;
    p.Nresidual = 1e1;
else
    p.Nbatch    = 1;
    p.Nperm     = 1e2;
    p.Nresidual = 1e3;
    p.permMethod       = 'FY_double_for_loop'; %p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
    p.permActionMethod = 'index';
end

names    = {'uniform','optimal','random'};
marksize = [15 8 8];
for kk=3:-2:1
    for pind=1:nout(kk)
        % construct measurement grid
        if kk==1
            [t,~] 
        elseif kk==2
            t = xmaster(pind,:);
        elseif kk==3
            t = rand(1,p.Nmeas);
        end
        hvbool='off'; % decide if it should be in legend
        if pind==1
            hvbool='on';
        end
        Jvec = wrap_getCostFun(t,p,active_inds); % get cost fun
%         Jvec = wrap_twoper_cfun(t,p,active_inds,f1,f2);
        % plot
        for ii=1:4
            for jj=1:ii
                nexttile(jj + 4*(ii-1))
                loglog(Jvec(jj),Jvec(ii+1),'.k','Color',cols(kk), ...
                    'MarkerSize',marksize(kk), ...
                    'DisplayName',names{kk},'HandleVisibility',hvbool);
                hold on
                xlabel(lab_list{jj})
                ylabel(lab_list{ii+1})
                drawnow
            end
        end
    end
end

%%

linkaxes([nexttile(5) nexttile(6)],'y')
linkaxes([nexttile(9) nexttile(10) nexttile(11)],'y')
linkaxes([nexttile(13) nexttile(14) nexttile(15) nexttile(16)],'y')

linkaxes([nexttile(1) nexttile(5) nexttile(9) nexttile(13)],'x')
linkaxes([nexttile(6) nexttile(10) nexttile(14)],'x')
linkaxes([nexttile(11) nexttile(15)],'x')
% %%
% for ii=1:4
%     for jj=1:ii
%         nexttile(jj + 4*(ii-1))
%         xl=xlim;
%         xticks([linspace(xl(1),0.5*xl(2),3)])
%         round(xticks,1)
%         tix=get(gca,'xtick')';
%         set(gca,'xticklabel',num2str(tix,'%0.1f'))
%     end
% end
% 
% for ii=1:4
%     for jj=1:ii
%         nexttile(jj + 4*(ii-1))
%         yl=ylim;
%         yticks(yl)
%         tix=get(gca,'ytick')';
%         set(gca,'yticklabel',num2str(tix,'%0.1f'))
%     end
% end

%%

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
leg=legend('Location','southoutside','NumColumns',3)
leg.Layout.Tile='south'
%%
set(findall(gcf,'-property','FontSize'),'FontSize',7)
print(gcf,strcat('/home/turner/research/overleaf/guelph_talk/figs_csc/',fname),'-dpng','-r600') % -r sets the resolution
savefig(fname)

function Jvec = wrap_twoper_cfun(t,p,active_inds,f1,f2)
p.freq=f1;
Jvec1=wrap_getCostFun(t,p,active_inds);
p.freq=f2;
Jvec2=wrap_getCostFun(t,p,active_inds);

Jvec=0.5*(Jvec1+Jvec2);
end