syms k lambda
solve(diff(exp(-lambda/2)*(lambda/2)^k/factorial(k),lambda),lambda)
%%
simplify(diff(diff(exp(-lambda/2)*(lambda/2)^k/factorial(k),lambda),lambda))

%%
close all
Sk = @(lambda,k)exp(-lambda/2).*(lambda/2).^k/factorial(k);

xv=0:.1:30;
plot(xv,Sk(xv,1))
xline(2)

%%
ncfsum = @(x,n1,n2,lambda,nn,NN) sum(exp(-lambda/2)*(0.5*lambda).^(nn:NN).*arrayfun(@(jj) regbinc(n1*x/(n2+n1*x),n1/2+jj,n2/2),nn:NN)./factorial(nn:NN));
regbinc=@(x,a,b) betainc(x,a,b)/betainc(1,a,b);

close all
xv=5:.1:50;
plot(xv,arrayfun(@(x) ncfsum(x,2,N-3,5,0,2),xv))
hold on
plot(xv,ncfcdf(xv,2,N-3,1))


lambda=1
%%
close all
M=1;
x0=rand
xv=0:.1:100;
plot(xv(1:end-1),diff(arrayfun(@(x) ncfsum(x0,2,N-3,x,0,M),xv))/.1)
hold on
plot(xv(1:end-1),-0.5*arrayfun(@(x) ncfsum(x0,2,N-3,x,M,M),xv(1:end-1)))

%% CHeck series convergence
N=8;
x=rand;
n1=2;
n2=N-3;
%%
gk=@(kk) exp(-kk)*kk^kk*regbinc(n1*x/(n2+n1*x),n1/2+kk,n2/2)/factorial(kk);
K=10
sum(arrayfun(@(kk)gk(kk),1:K))
%%
u=8
.1^u*u^3