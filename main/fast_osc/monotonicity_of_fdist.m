N=8;
n1=2;
n2=N-3;
%% check monotonicity naively
N=randsample(4:15,1);
xv=0:.001:1;
x=100*rand;
prod(diff(ncfcdf(x,2,N-3,xv))<0)
%% define beta function and terms Tk 
regbinc=@(x,a,b) betainc(x,a,b)/betainc(1,a,b);

ncfsum = @(x,n1,n2,lambda,nn,NN) sum(exp(-lambda/2)*(0.5*lambda).^(nn:NN).*arrayfun(@(jj) regbinc(n1*x/(n2+n1*x),n1/2+jj,n2/2),nn:NN)./factorial(nn:NN));
T      = @(x,n1,n2,lambda,k) exp(-lambda/2).*(0.5*lambda).^k*regbinc(n1*x/(n2+n1*x),n1/2+k,n2/2)/factorial(k);
g      = @(x,n1,n2,k) exp(-k).*(k.^k)*regbinc(n1*x/(n2+n1*x),n1/2+k,n2/2)/factorial(k);

close all

%% check I{k+1} <= Ik 
close all
N=randsample(4:14,1);
n1=2;
n2=N-3;

x=0:.1:30;
for k=1:10
    plot(x,regbinc(n1*x./(n2+n1*x),n1/2+k,n2/2))
    hold on
    pause(0.1)
end
%% ratio test for convergence in lambda
N=randsample(4:14,1);

lambda=10*rand;       % non-centrality
k=randsample(0:20,1); % randomly choose one of first twenty terms

% need to assume x > 0 or else have issue at end point
xv=0:.01:10;
Tkp1 = arrayfun(@(x) T(x,2,N-3,lambda,k+1),xv);
Tk   = arrayfun(@(x) T(x,2,N-3,lambda,k),xv);
rat  =Tkp1./Tk ;

% check ratio for all x except end point
prod(rat(2:end) < lambda/(k+1)/2)

%% check max occurs at lambda=2k
syms k lambda
% critical points
solve(diff(exp(-lambda/2)*lambda^k/factorial(k),lambda)==0) 
% check second derivative is negative at lambda=2k 
simplify(subs(diff(diff(exp(-lambda/2)*lambda^k/factorial(k),lambda),lambda),lambda,2*k))

%% termwise derivative bound
lambda=10*rand;       % non-centrality
k=randsample(0:20,1); % randomly choose one of first twenty terms
h=.001;
lv=0.1:h:50;

dTk = diff(T(x,2,N-3,lv,k))/h
min(abs(dTk))
% fails when lambda << 1 but probably to be expected since |dTk| ~ 1e-28
prod(dTk<-0.5*T(x,2,N-3,lv(1:end-1),k) + 0.5*T(x,2,N-3,lv(1:end-1),k-1))

%% global derivative bound
clf
x=10*rand;
k=randsample(4:20,1);
h=.001
lv=0:h:50

dTk = diff(T(x,2,N-3,lv,k))/h;
plot(lv(1:end-1),abs(dTk))
hold on
Gk  = 0.5*abs(g(x,2,N-3,k) + g(x,2,N-3,k-1) );
yline(Gk)

%% verify ratio of g_{k+1} g_k is correct
n1=2;n2=N-3;
x = 10*rand;
k = randsample(4:10,1);

Ik = @(k) regbinc(n1*x/(n2+n1*x),n1/2+k,n2/2);
fprintf("%.8f\n",g(x,2,N-3,k+1)/g(x,2,N-3,k))
fprintf("%.8f\n",(Ik(k+1)/Ik(k))*(1+1/k)^k*exp(-1))


%% visually check uniform convergence
close all
N=randsample(4:15,1);
lv=0:.1:40;
x=rand*10
F=ncfcdf(x,n1,n2,lv);
plot(lv,ncfcdf(x,2,N-3,lv),'-k') % true cdf
hold on

for ii=1:10
    nstop=ii;
    clf
    S=arrayfun(@(ll) ncfsum(x,2,N-3,ll,0,nstop),lv);
    plot(lv,ncfcdf(x,2,N-3,lv),'-k') % true cdf
    hold on
    plot(lv,S,'--k')
    xlabel('lambda')
    pause(0.1)
end
% check that derivative is dominated by zeroth term of function
% close all
% x=10*rand
% T0= @(x,n1,n2,lambda) -1/2*exp(-lambda/2)*regbinc(n1*x./(n2+n1*x),n1/2,n2/2)
% 
% h=.01;
% lv=0:h:100;
% F=ncfcdf(x,n1,n2,lv);
% Fp=diff(F)/h;
% 
% plot(lv(1:end-1),Fp,'-k')
% hold on
% plot(lv,T0(x,n1,n2,lv),'--k')


%%
close all
plot(1,1)
%%
h=.1;
x=rand*20;
N=8;
close all
for k=1:20
lv=0:h:100;
F=ncfcdf(x,2,N-3,lv);
Fp=diff(F)/h;

clf
Sk = arrayfun(@(ll) ncfsum(x,2,N-3,ll,0,k),lv);
Tk = arrayfun(@(ll) ncfsum(x,2,N-3,ll,k,k),lv);
Skp = diff(Sk)/h;

plot(lv(1:end-1),Fp,'-r')
hold on
plot(lv(1:end-1),Skp,'-k')

plot(lv,-0.5*Tk,'--k')
ylim([-.1,0])
pause(0.08)
end
%% Check for small NN
close all
T0 = arrayfun(@(ll) ncfsum(x,2,N-3,ll,0,0),0:h:20);
T1 = arrayfun(@(ll) ncfsum(x,2,N-3,ll,1,1),0:h:20);
T2 = arrayfun(@(ll) ncfsum(x,2,N-3,ll,2,2),0:h:20);

T0p = diff(T0)/h;
T1p = diff(T1)/h;
T2p = diff(T2)/h;

plot(T0p+T1p+T2p,'-k')
hold on
plot(-0.5*T2,'--k')
%% Check sums satisfy bound
close all
k=2;


