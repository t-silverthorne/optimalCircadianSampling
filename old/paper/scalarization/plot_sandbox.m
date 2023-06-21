load('run_sandbox_out_amp_higher.mat')
figure
addpath('../utils_core/')
addpath('../utils_cost_fun/')
tiledlayout(3,2)

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');    
set(groot,'defaultLegendInterpreter','latex');

xlab_list={'amp bias','acro bias','amp variance', 'acro variance', '1-power'};


[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
J_unif = wrap_getCostFun(t_unif,p,active_inds);


disp("finished1")

for ii=1:size(xmaster,1)
%     t=NaN(1,length(xmaster(ii,:))+1);
%     t(1)=0;
%     for jj=2:length(t)
%         t(jj) = t(jj-1) + xmaster(ii,jj-1)*(1-t(jj-1));
%     end
    t=xmaster(ii,:);
    J=wrap_getCostFun(t,p,active_inds);
    
    for jj=1:5
        nexttile(jj)    
        if jj==5
            plot(J(jj),J(1),'.k','MarkerSize',7);
            hold on
            if (ii==1 || ii == size(xmaster,1))
               plot(J_unif(jj),J_unif(1),'.r','MarkerSize',10);
            end
        else
            plot(J(jj),J(jj+1),'.k','MarkerSize',7);
            hold on
            if (ii==1 || ii == size(xmaster,1))
               plot(J_unif(jj),J_unif(jj+1),'.r','MarkerSize',20);
            end
        end
        
    end

    for jj=1:5
        nexttile(jj)
        xlabel(xlab_list{jj})
        if jj==5
            ylabel(xlab_list{1});
        else
            ylabel(xlab_list{jj+1});
        end
    end
    drawnow
   
end
disp("finished2")
