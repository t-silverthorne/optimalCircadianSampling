% path and data
load('data/test_amp4.mat')
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

% tiling
tiledlayout(2,3,'TileSpacing','tight')

% actual code
nout(1)=1;
nout(2)=size(xmaster,1);
nout(3)=size(xmaster,1);

nin = 1;
c1="#648FFF"; % uniform colour
c2="#785EF0"; % 2-uniform colour
c3="#DC267F"; % Unif(0,1) colour
c4="#FE6100"; % jittered colour
c5="#000000";
cols=[c1,c5,c3];

p.Nperm=1e2;
p.Nresidual=1e1;
ntypes=3;
names={'uniform','optimal','random'};

for jj=1:ntypes

    for pind=1:nout(jj)
        if jj==1
            [t,~] = getSamplingSchedules(p.Nmeas,0,0,0);
        elseif jj==2
            t = xmaster(pind,:);
        elseif jj==3
            t = rand(1,p.Nmeas);
        end
        hvbool='off';
        if pind==1
            hvbool='on';
        end
        Jvec = wrap_getCostFun(t,p,active_inds);
    
        nexttile(1)
        plot(Jvec(1),Jvec(3),'.k','Color',cols(jj),'MarkerSize',10, ...
            'DisplayName',names{jj},"HandleVisibility",hvbool)
        hold on
        nexttile(2)
        plot(Jvec(4),Jvec(3),'.k','Color',cols(jj),'MarkerSize',10, ...
            'DisplayName',names{jj},"HandleVisibility",hvbool)
        hold on
        nexttile(3)
        plot(Jvec(5),Jvec(3),'.k','Color',cols(jj),'MarkerSize',10, ...
            'DisplayName',names{jj},"HandleVisibility",hvbool)
        hold on
        
        nexttile(4)
        plot(Jvec(1),Jvec(4),'.k','Color',cols(jj),'MarkerSize',10, ...
            'DisplayName',names{jj},"HandleVisibility",hvbool)
        hold on
        nexttile(5)
        plot(Jvec(2),Jvec(4),'.k','Color',cols(jj),'MarkerSize',10, ...
            'DisplayName',names{jj},"HandleVisibility",hvbool)
        hold on
        nexttile(6)
        plot(Jvec(5),Jvec(4),'.k','Color',cols(jj),'MarkerSize',10, ...
            'DisplayName',names{jj},"HandleVisibility",hvbool)
        hold on

        drawnow
    end
end

% more aesthetics 
nexttile(1)
xlabel('amp bias')
ylabel('amp var')
set(gca,'XTickLabel',[]);

nexttile(2)
xlabel('1+acro var')
set(gca,'YTickLabel',[]);

nexttile(3)
xlabel('power')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


nexttile(4)
xlabel('amp bias')
ylabel('acro var')

nexttile(5)
xlabel('acro bias')
set(gca,'YTickLabel',[]);

nexttile(6)
xlabel('power')
set(gca,'YTickLabel',[]);

nexttile(5)
legend('Location','southoutside','NumColumns',3)