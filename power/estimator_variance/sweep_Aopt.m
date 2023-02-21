clear
close all
addpath('../utils/')
param.useGPU=false;
param.NL=8;
param.NR=5;
param.useTraceMC=false;
param.NumTraceMC=1e4;

xvar='freq_true';
cvar='NL';
xmin=1;
xmax=8;
Numx=1e3;

xvals=linspace(xmin,xmax,Numx);
cmin=1;
cmax=10;
cvals=cmin:cmax;
Numc=numel(cvals);

blue=[0 0.44 1];
red =[0.8 0 0];
color_val=@(ind,Numc) red*(ind-Numc)/(1-Numc) + blue*(1-ind)/(1-Numc);

for ii=1:Numc
    yvals=NaN(1,Numx);
    param.(cvar)=cvals(ii);
    [~,t]=getSamplingSchedules(param.NL,param.NR,0,0.5);
    for jj=1:Numx    
        param.(xvar)=xvals(jj);
        yvals(jj)=get_Aoptimal_cost(t,param);
    end
    semilogy(xvals,yvals,'-k','color',color_val(ii,Numc))
    hold on
    pause(0.5)
end

%% v2
close all
clf

xvar='freq_true';
cvar='tau';
pvar='NL';
xmin=1;
xmax=12;
Numx=1e3;
xvals=linspace(xmin,xmax,Numx);

cmin=0.5;
cmax=0.25;
Numc=10;
cvals=linspace(cmin,cmax,Numc);

pmin=4;
pmax=8;
pvals=pmin:1:pmax;
Nump=numel(pvals);
tiledlayout(1,Nump)

for pp=1:Nump
    nexttile(pp)
    param.NL=pvals(pp);
    param.NR=pvals(pp);
    for ii=1:Numc
        yvals=NaN(1,Numx);
        tau=cvals(ii);
        [~,t]=getSamplingSchedules(param.NL,param.NR,0,tau);
        for jj=1:Numx    
            param.(xvar)=xvals(jj);
            yvals(jj)=get_Aoptimal_cost(t,param);
        end
        semilogy(xvals,yvals,'-k','color',color_val(ii,Numc))
        hold on
    end
end


function cost = get_Aoptimal_cost(t,param)
X=constructX(t,param);
if param.useTraceMC
    d=size(X'*X,1);
    eps=randn(d,1,param.NumTraceMC);
    cost=sum(pagemtimes(pagetranspose(eps),pagemldivide(X'*X,eps)),3)/param.NumTraceMC;
else
    cost = trace(inv(X'*X));
end
end