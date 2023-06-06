clf
clear
load('data/test_amp4.mat')
fname='test_amp4fig';
testing=true;

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


lab_list = {'amp bias','acro bias','amp var','1-acro var', '1-power'};

tiledlayout(4,4,'TileSpacing','none')
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

if testing
    p.Nbatch    = 1;
    p.Nperm     = 1e1;
    p.Nresidual = 1e1;
else
    p.Nbatch    = 1;
    p.Nperm     = 1e2;
    p.Nresidual = 1e3;
end

names    = {'uniform','optimal','random'};
marksize = [20 10 10];
for kk=3:-1:1
    for pind=1:nout(kk)
        % construct measurement grid
        if kk==1
            [t,~] = getSamplingSchedules(p.Nmeas,0,0,0);
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



linkaxes([nexttile(5) nexttile(6)],'y')
linkaxes([nexttile(9) nexttile(10) nexttile(11)],'y')
linkaxes([nexttile(13) nexttile(14) nexttile(15) nexttile(16)],'y')

linkaxes([nexttile(1) nexttile(5) nexttile(9) nexttile(13)],'x')
linkaxes([nexttile(6) nexttile(10) nexttile(14)],'x')
linkaxes([nexttile(11) nexttile(15)],'x')


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