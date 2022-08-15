% compare Wynn's method and modified Wynn's method when residuals have
% memory

% tv1=sort(rand(d*2,1));
% tv2=sort(rand(d,1));
% wv1=ones(d*2,1)/numel(tv1);
% wv2=ones(d,1)/numel(tv2);
% 
% get_H2_matrix(tv2,wv2,tv1,wv1,3,'exp')


clear
clf
sigval=.1;
gamval=10;
niter=500;
set_memory_params(sigval,gamval)


orange=[0.8500 0.3250 0.0980];
blue=[0 0.4470 0.7410];


%rng(1) % nice example of convergence
k=5;
kernel_opt='exp';

cvals=0;
for i=1:niter
    theta_vec=linspace(0,1,i+1);
    theta_vec=theta_vec(1:end-1)';
    w_vec=ones(numel(theta_vec),1)/numel(theta_vec);
    cvals(i)=log(det(get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt)));
end
%semilogy(1:niter,cvals,'.k')


% theta_vec=linspace(0,1,5+1);
% theta_vec=theta_vec(1:end-1)';
% w_vec=ones(numel(theta_vec),1)/numel(theta_vec);
% log(det(get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt)))
 
[~,mind]=max(real(cvals))
display(max(real(cvals)))
%%

rng(2) % used in presentation
n0=1;
scur=rng;
clf
theta_vec=rand(n0,1);
w_vec=ones(n0,1)/numel(theta_vec);

warning('off','MATLAB:nearlySingularMatrix')


h=NaN(1,3);
for i=1:niter
    Mnow=get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt);
    f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];
    H1=@(x) get_H1_matrix(theta_vec,w_vec,x,1,k,kernel_opt);
    H2=@(x) get_H2_matrix(theta_vec,w_vec,x,1,k,kernel_opt);
    xnew=fminbnd(@(x) -trace(Mnow\H1(x) -2*Mnow\H2(x)),0,1);
    theta_vec=[theta_vec; xnew];
    w_vec=ones(n0+i,1)/numel(theta_vec);
    semilogy(i,log(det(get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt))), ...
    '.k','MarkerEdgeColor',orange')
    if i==niter
        h(1)=semilogy(i,log(det(get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt))), ...
    '.k','MarkerEdgeColor',orange','DisplayName','memory Wynn');
    end
    hold on
end

n0=1;
scur=rng;
theta_vec=rand(n0,1);
w_vec=ones(n0,1)/numel(theta_vec);

for i=1:niter
    Mnow=get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt);
    f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];
    xnew=fminbnd(@(x) -trace(Mnow\get_memory_info_matrix_k(x,1,k,kernel_opt)),0,1);
    theta_vec=[theta_vec; xnew];
    w_vec=ones(n0+i,1)/numel(theta_vec);
    semilogy(i,log(det(get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt))), ...
    'xk','MarkerEdgeColor','black','MarkerSize',1.5)
    if i==niter
        h(2)=semilogy(i,log(det(get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt))), ...
    'xk','MarkerEdgeColor','black','DisplayName','original Wynn','MarkerSize',1.5);
    end
    hold on
end


warning('on','MATLAB:nearlySingularMatrix')
%ylim([-10 0])
h(3)=yline(max(real(cvals)),'--k','color',blue,'DisplayName','optimal uniform','LineWidth',1.3)
legend(h,'location','southeast')
xlabel('iteration')
ylabel('$\log(\det(M(\xi_{n}))$','Interpreter','latex')
xticks([0 100 200 300 400 500])

%%

ylabel('$\Psi(\xi_{n})$','Interpreter','latex')
ylim([10^0 25])
yticks([10^0 10^1 2*10^1])
yticklabels({'1','$10$','$20$'})

plot_filename='modified_wynn'
ht=2; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too


%% Compare to simulated annealing

%%
N=30; % number of sampling points
w=rand(1,N);
w=w/sum(w);
options=optimoptions('simulannealbnd','PlotFcns',{'saplotf','saplotbestf'}, ...
                     'MaxIterations',2000);
[xval,fval]=simulannealbnd( @(x) sa_info_matrix_cost_fun([x(1:N) x((N+1):end)/sum(x((N+1):end))],k), ...
    rand(1,2*N),zeros(1,2*N),[ones(1,N) ones(1,N)],options)


function C = sa_info_matrix_cost_fun(x,k)
memory=true;
n_over2=numel(x)/2;
theta_vec=x(1:n_over2)';
w_vec=x((n_over2+1):end)';
if memory % take real part, necessary if using degenerate number of sampling points
    C=real(-log(det(get_memory_info_matrix_k(theta_vec,w_vec,k,'exp'))));
else
    C=-log(det(get_info_matrix(theta_vec,w_vec)));
end
end


function M = get_memory_info_matrix_k(theta_vec,w_vec,k,kernel_opt)
f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];
global sigma0
global gamma0

switch kernel_opt
    case 'exp'
        alpha = @(x,xp) exp(-abs(.15*(x-xp)));
    case 'tanhsq'
        alpha = @(x,xp) 1-tanh(abs(.15*(x-xp))).^2;
end


M=0;
for i=1:numel(theta_vec)
    ss=sigma0*(1+ gamma0*sum(w_vec(1:i-1).*alpha(theta_vec(i),theta_vec(1:i-1))));
    M=M+w_vec(i)*f(theta_vec(i))'*f(theta_vec(i))/ss^2;
end
end


function H1 = get_H1_matrix(theta_vec_star,w_vec_star,theta_vec,w_vec,k,kernel_opt)
f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];
global sigma0
global gamma0

switch kernel_opt
    case 'exp'
        alpha = @(x,xp) exp(-abs(.15*(x-xp)));
    case 'tanhsq'
        alpha = @(x,xp) 1-tanh(abs(.15*(x-xp))).^2;
end


H1=0;
for ii=1:numel(theta_vec)
    jj=find(theta_vec_star<theta_vec(ii),1,'last');
    ss=sigma0*(1+ gamma0*sum(w_vec_star(1:jj).*alpha(theta_vec_star(jj),theta_vec_star(1:jj))));
    H1=H1+w_vec(ii)*f(theta_vec(ii))'*f(theta_vec(ii))/ss^2;
end

end

function H2 = get_H2_matrix(theta_vec_star,w_vec_star,theta_vec,w_vec,k,kernel_opt)
f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];

global sigma0
global gamma0

switch kernel_opt
    case 'exp'
        alpha = @(x,xp) exp(-abs(.15*(x-xp)));
    case 'tanhsq'
        alpha = @(x,xp) 1-tanh(abs(.15*(x-xp))).^2;
end

H2=0;
for ii=1:numel(theta_vec)
    jj=find(theta_vec_star<theta_vec(ii),1,'last');
    ss_star=sigma0*(1+ gamma0*sum(w_vec_star(1:jj).*alpha(theta_vec_star(jj),theta_vec_star(1:jj))));
    ss_minus_1=sigma0*(gamma0*sum(w_vec(1:ii-1).*alpha(theta_vec(ii),theta_vec(1:ii-1))));
    
    H2=H2+w_vec(ii)*f(theta_vec(ii))'*f(theta_vec(ii))*ss_minus_1/ss_star^3;
end
end

function set_memory_params(sigval,gamval)
global sigma0
global gamma0
sigma0=sigval;
gamma0=gamval;
end