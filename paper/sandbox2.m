close all
p.Nacro     = 1;
p.Nresidual = 1;
p.freq      = 7.8;
p.Amp       = 1.5;
p.noise     = .8;
p.Nmeas     = 16;

t_hires=linspace(0,1,100);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);

Y=getSimulatedData(t_unif,p);
%Y=Y(:,:,:,acind);
[pxx1,f1]=periodogram(Y',[],[],p.Nmeas);
[pxx2,f2]=periodogram(Y',t_unif);
[pxx3,f3]=plomb(Y',t_unif,[],30);
plot(f1,pxx1,'.-r')
hold on
%plot(f2*p.Nmeas/pi/2,pxx2,'--k')
plot(f3,pxx3,'.-k')