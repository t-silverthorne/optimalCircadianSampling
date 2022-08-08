% 
clear
clf
orange=[0.8500 0.3250 0.0980];
blue=[0 0.4470 0.7410];

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
T1est=23/24;
T2est=5/24;
tv=get_theta_vec_atwood(5,100,T1est,T2est);
theta_unif=linspace(0,1)';
theta_unif=theta_unif(1:end-1);
w_vec_unif=ones(numel(theta_unif),1)/numel(theta_unif);
log(det(get_info_matrix(theta_unif,w_vec_unif,T1est,T2est)))
w_vec=ones(numel(tv),1)/numel(tv);
log(det(get_info_matrix(tv,w_vec,T1est,T2est)))

function theta_vec= get_theta_vec_atwood(n0,niter,T1,T2)

theta_vec=rand(n0,1);
w_vec=ones(n0,1)/numel(theta_vec);

gamma=@(n) (n+1)^(-1);

tiledlayout(2,1)

for i=1:niter
    Mnow=get_info_matrix(theta_vec,w_vec,T1,T2);

    f=@(theta) [1 cos(2*pi*theta/T1) sin(2*pi*theta/T1) cos(2*pi*theta/T2) sin(2*pi*theta/T2)];
    [xnew_best,fnew_plus]=fminbnd(@(x) -trace(Mnow\(f(x)'*f(x))),0,1); % get best point to add
    [fold_worst ,ind_bad]=max( arrayfun(@(ind) ...
        trace(Mnow\(f(theta_vec(ind)))'*f(theta_vec(ind))), ...
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
end
end


function M = get_info_matrix(theta_vec,w_vec,T1,T2)
x=@(theta) [1 cos(2*pi*theta/T1) sin(2*pi*theta/T1) cos(2*pi*theta/T2) sin(2*pi*theta/T2)];

M=0;
for i=1:numel(theta_vec)
    M=M+w_vec(i)*x(theta_vec(i))'*x(theta_vec(i));
end
end