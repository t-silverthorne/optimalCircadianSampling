function [acrovec,pwr]=reshaped_simulatePWR_fully_vectorized(param,nodes)
NL=param.NL;
NR=param.NR;
Nperm=param.Nperm;
Nresidual=param.Nresidual;
Nacro=param.Nacro;
Amp=param.Amp;
Nmeas=NL+NR;
freq_true=param.freq_true; % do not confuse with freq_est

if isnumeric(nodes)
    t=nodes;
    Nmeas=length(t);
else
    switch nodes
        case 'uniform'
            [t,~]=getSamplingSchedules(NL,NR,0,0.5);
        case 'non-uniform'
            [~,t]=getSamplingSchedules(NL,NR,0,0.3);
    end

end
X=constructX(t,param); % construct linear model
acrovec=linspace(0,2*pi,Nacro);
acromat=reshape(acrovec,1,1,1,Nacro);

eps=randn(Nmeas,Nresidual,1,Nacro,'gpuArray');
use_same_perms=true;

if use_same_perms
    permMat=rand(Nmeas,Nresidual,Nperm,'gpuArray');
    [~,I]=sort(permMat,1);
    I=repmat(I,1,1,1,Nacro);
else
    permMat=rand(Nresidual,Nmeas,Nperm,Nacro,'gpuArray');
    [~,I]=sort(permMat,2);
end

t=transpose(t);
clear permMat
Y=Amp*cos(2*pi*t*freq_true-acromat)+param.noise*eps;
clear eps


betas_obs=pagemldivide(X'*X,pagemtimes(X',Y)); % observed error

fits_obs=pagemtimes(X,betas_obs);
SSres_obs=sqrt(sum((fits_obs-Y).^2,2));
clear fits_obs

n=size(Y,1);m=size(Y,2);r=size(Y,3);s=size(Y,4);
inds=0:m*n*r:m*n*r*s;
inds=inds(1:end-1);
inds=reshape(inds,1,1,1,length(inds));
offMat=repmat((0:m-1)',1,n)*n;

YI=Y(pagetranspose(I+offMat+inds));
clear I offMat inds
%max(max(max(max(abs(sort(YI,2)-sort(Y,2))))))

betas=pagemldivide(X'*X,pagemtimes(X',YI));

fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr=sum(sum(SSres>SSres_obs,3)/Nperm>.95)/Nresidual;
pwr=reshape(pwr,1,Nacro);
end