% ignore problem structure, just do naive approach to minimize -det(M)
% using simulated annealing
addpath('../utils')
clear
N=3; % number of sampling points
w=rand(1,N);
w=w/sum(w);
options=optimoptions('simulannealbnd','PlotFcns',{'saplotf','saplotbestf'}, ...
                     'MaxIterations',1000);
simulannealbnd( @(x) sa_info_matrix_cost_fun([x(1:N) x((N+1):end)/sum(x((N+1):end))]), ...
    rand(1,2*N),zeros(1,2*N),[ones(1,N) ones(1,N)],options)
%%

theta_unif=linspace(0,1,6+1)';
theta_unif=theta_unif(1:end-1);
w_vec_unif=ones(numel(theta_unif),1)/numel(theta_unif);
-log(det(get_memory_info_matrix(theta_unif,w_vec_unif,'exp')))


function C = sa_info_matrix_cost_fun(x)
memory=true;
n_over2=numel(x)/2;
theta_vec=x(1:n_over2);
w_vec=x((n_over2+1):end);
if memory
    C=-log(det(get_memory_info_matrix(theta_vec,w_vec,'exp')));
else
    C=-log(det(get_info_matrix(theta_vec,w_vec)));
end
end
