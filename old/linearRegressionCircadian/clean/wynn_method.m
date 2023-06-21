clear
orange=[0.8500 0.3250 0.0980];
blue=[0 0.4470 0.7410];

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Illustrates non-monotonic convergence
rng(1) % nice example of convergence
k=6;
n0=1;
theta_vec=rand(n0,1);
w_vec=ones(n0,1)/numel(theta_vec);

h=NaN(1,3);
for i=1:100
    Mnow=get_info_matrix_k(theta_vec,w_vec,k);
    f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];
    xnew=fminbnd(@(x) -trace(Mnow\(f(x)'*f(x))),0,1);
    theta_vec=[theta_vec; xnew];
    w_vec=ones(n0+i,1)/numel(theta_vec);
    semilogy(i,log(det(get_info_matrix_k(theta_vec,w_vec,k))),'.k','MarkerFaceColor',orange,'MarkerEdgeColor',orange)
    if i==100
        h(1)=semilogy(i,log(det(get_info_matrix_k(theta_vec,w_vec,k))),'.k', ...
            'MarkerFaceColor',orange,'MarkerEdgeColor',orange, 'DisplayName','$\log(\det(M(\xi_{n})))$');
    end
    hold on
end
h(3)=xline(k*(k+1)/2,'DisplayName','$k(k+1)/2$');


theta_unif=linspace(0,1,8*k)';
theta_unif=theta_unif(1:end-1);
w_vec_unif=ones(numel(theta_unif),1)/numel(theta_unif);
h(2)=yline(log(det(get_info_matrix_k(theta_unif,w_vec_unif,k))),'--k', ...
    'color',blue,'DisplayName','$\log(\det(M(\xi^{*})))$');

hold off
%%
legend(h,'location','southeast')
ylabel('$\log(\det(M(\xi_{n})))$')
xlabel('$n$')

plot_filename='wynn_convergence';
ht=2.5; % height
wd=6; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too

function M = get_info_matrix_k(theta_vec,w_vec,k)
syms theta
assume(theta,'real')
f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];

M=0;
for i=1:numel(theta_vec)
    M=M+w_vec(i)*f(theta_vec(i))'*f(theta_vec(i));
end
end
