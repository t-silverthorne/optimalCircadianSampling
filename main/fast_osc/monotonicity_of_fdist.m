%% check regularized beta follow proper pattern
regbinc=@(x,a,b) betainc(x,a,b)/betainc(1,a,b);

%%
close all
N=8;
n1=2;
n2=N-3;
x=0:.1:30;
for k=1:10
    plot(x,regbinc(n1*x./(n2+n1*x),n1/2+k,n2/2))
    hold on
    pause(0.1)
end

%% check that derivative is dominated by zeroth term of function
close all
x=10*rand
T0= @(x,n1,n2,lambda) -1/2*exp(-lambda/2)*regbinc(n1*x./(n2+n1*x),n1/2,n2/2)

h=.01;
lv=0:h:100;
F=ncfcdf(x,n1,n2,lv);
Fp=diff(F)/h;

plot(lv(1:end-1),Fp,'-k')
hold on
plot(lv,T0(x,n1,n2,lv),'--k')

%% check monotonic naively
N=randsample(4:15,1);
xv=0:.001:1;
x=100*rand;
prod(diff(ncfcdf(x,2,N-3,xv))<0)

%% check term-wise derivatives satisfy bound

ncfsum = @(x,n1,n2,lambda,nn,NN) sum(exp(-lambda/2)*(0.5*lambda).^(nn:NN).*arrayfun(@(jj) regbinc(n1*x/(n2+n1*x),n1/2+jj,n2/2),nn:NN)./factorial(nn:NN));
close all
x=10*rand
N=8;

h=.1;
k=3;
Tk=arrayfun(@(ll) ncfsum(x,2,N-3,ll,k,k),0:h:20);
Tkm1=arrayfun(@(ll) ncfsum(x,2,N-3,ll,k-1,k-1),0:h:20);

plot(diff(Tk)/h,'-k')
hold on
plot( -0.5*Tk +0.5*Tkm1,'--k')



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
%%
close all
k=2;
Sk = arrayfun(@(ll) ncfsum(x,2,N-3,ll,0,k),0:h:20);
Tk = arrayfun(@(ll) ncfsum(x,2,N-3,ll,k,k),0:h:20);
Skp = diff(Sk)/h;

plot(Skp,'-k')
hold on
plot(-0.5*Tk,'--k')

%%
syms k lambda
diff(exp(-lambda/2)*lambda^k/factorial(k),lambda)