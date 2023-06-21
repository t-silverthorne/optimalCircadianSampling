% generates samples from pi(f), a truncated gaussian distribution
clear

close all
N=10^6;
muf=.1;
sigf=.1;
Zf0=(1125899906842624*2^(1/2)*pi^(1/2)*(erf((2^(1/2)*muf*(1/sigf^2)^(1/2))/2) + 1)*(sigf^2)^(1/2))/(5644425081792261*sigf);%
ftr = @(x) (x>0).*exp(-(x-muf).^2/2/sigf^2)/sigf/sqrt(2*pi)/Zf0;

samp=sampleTruncCDF(N,muf,sigf);
histogram(samp,floor(sqrt(length(samp))),'Normalization','pdf','EdgeColor','none')
xv=-1:.01:10;
hold on
plot(xv,ftr(xv))
xlim([0 1])


function x=sampleTruncCDF(N,muf,sigf)
Zf0=(1125899906842624*2^(1/2)*pi^(1/2)*(erf((2^(1/2)*muf*(1/sigf^2)^(1/2))/2) + 1)*(sigf^2)^(1/2))/(5644425081792261*sigf);%
Ftr = @(x) (x>0).*(erf((x-muf)/sigf/sqrt(2))  - erf((-muf)/sigf/sqrt(2)) )/2/Zf0;
ftr = @(x) (x>0).*exp(-(x-muf).^2/2/sigf^2)/sigf/sqrt(2*pi)/Zf0;

y=rand(N,1);

x=muf;%xv(ind);

for i=1:10
 x= x - (Ftr(x)-y)./ftr(x);
end
end