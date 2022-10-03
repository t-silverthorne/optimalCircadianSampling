addpath('utils/')

%% simple test
Ntimes=5
mt_unif=rand(1,Ntimes+1)*24;
mt_unif=mt_unif(1:end-1);

theta_vec=mt_unif/24;
w_vec=ones(Ntimes,1)/Ntimes;
M=0;
x=@(theta) [ones(numel(theta),1) cos(2*pi*theta) sin(2*pi*theta)];
for i=1:numel(theta_vec)
    M=M+w_vec(i)*x(theta_vec(i))'*x(theta_vec(i));
end
M
x(theta_vec')'*diag(w_vec)*x(theta_vec')

%% uniform example
x=@(theta) [ones(numel(theta),1) cos(2*pi*theta) sin(2*pi*theta)];

Nleft=10;
Nright=4;
method='lin';% nonlin or lin
Ntimes=Nleft+Nright;

mt1=linspace(0,8,Nleft+1);
mt1=mt1(1:end-1);
mt2=linspace(8,24,Nright+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2];

% construct uniform grid 
mt_unif=linspace(0,24,Ntimes+1);
mt_unif=mt_unif(1:end-1);

mt_nu=mt_nu/24;
mt_unif=mt_unif/24;
w_vec=ones(Ntimes,1)/Ntimes;

det(x(mt_unif')'*diag(w_vec)*x(mt_unif'))
det(x(mt_nu')'*diag(w_vec)*x(mt_nu'))

%% non-uniform example