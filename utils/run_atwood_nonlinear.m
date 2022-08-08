function [theta_vec,Psi_w_unif]=run_atwood_nonlinear(param,niter,n0,testing)
param.T1=param.T1/24;
param.T2=param.T2/24;
orange=[0.8500 0.3250 0.0980];
blue=[0 0.4470 0.7410];

warning('off','MATLAB:nearlySingularMatrix')
if nargin<2
    niter=100;
end
if nargin<3
    n0=1;
end
if nargin<4
    testing=false;
end

theta_vec=rand(n0,1);
w_vec=ones(n0,1)/numel(theta_vec);
gamma=@(n) (n+1)^(-1);
tiledlayout(2,1)

for i=1:niter
    Mnow=get_FIM_biharmonic_nonlin(theta_vec,w_vec,param);
    f_nonlinear=@(t,beta1,beta2,beta3,beta4,T1,T2)[1.0;sin((t.*pi.*2.0)./T1);cos((t.*pi.*2.0)./T1);sin((t.*pi.*2.0)./T2);cos((t.*pi.*2.0)./T2);1.0./T1.^2.*t.*pi.*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*2.0;1.0./T2.^2.*t.*pi.*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*2.0];
    f=@(theta) f_nonlinear(theta,param.beta1,param.beta2,param.beta3,param.beta4,param.T1,param.T2);
    [xnew_best,fnew_plus]=fminbnd(@(x) -trace(Mnow\(f(x)*f(x)')),0,1); % get best point to add
    [fold_worst ,ind_bad]=max( arrayfun(@(ind) ...
        trace(Mnow\(f(theta_vec(ind)))*f(theta_vec(ind))'), ...
        1:numel(theta_vec))); % get worst point you have
    if fnew_plus<-fold_worst
        xnew=xnew_best;
        theta_vec=[theta_vec; xnew];
        w_vec=[(1-gamma(n0+i))*w_vec; gamma(n0+i)];
    else
        as=-min(gamma(n0+i),w_vec(ind_bad)/(1-w_vec(ind_bad)));
        w_vec=(1-as)*w_vec;
        w_vec(ind_bad)=w_vec(ind_bad)+as;
    end

    if testing
        
        if i==niter
            h(1)=plot(n0+i,log(det(get_FIM_biharmonic_nonlin(theta_vec,w_vec,param))), ...
                '.k','DisplayName','Atwood method','MarkerSize',5);
        else
            plot(n0+i,log(det(get_FIM_biharmonic_nonlin(theta_vec,w_vec,param))), ...
            '.k','MarkerSize',5)
        end
        hold on
    end
end
Psi_w_unif=log(det(get_FIM_biharmonic_nonlin(theta_vec, ...
    ones(numel(theta_vec),1)/numel(theta_vec),param)));

if testing
    theta_unif=linspace(0,1,numel(theta_vec)+1)';
    theta_unif=theta_unif(1:end-1);
    w_vec_unif=ones(numel(theta_unif),1)/numel(theta_unif);
    c_unif=log(det(get_FIM_biharmonic_nonlin(theta_unif,w_vec_unif,param)));
    h(2)=yline(c_unif,'--k','color',blue,'DisplayName','uniform design');
    hold off
    xlabel('iteration')
    ylabel('$\log(\det(M(\xi_{n}))$','Interpreter','latex')
    ylim([-100 20])
    legend(h,'Location','southeast','FontSize',12)


    plot_filename='atwood_design'
    ht=2; % height
    wd=6; % width
    set(gcf,'PaperUnits','inches')
    set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
    print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
    savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too


end

warning('on','MATLAB:nearlySingularMatrix')
end
