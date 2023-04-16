tic
p.Nmeas=10;
p.Nacro=16;
p.Nresidual=1e3;
p.Nperm=1e2;
p.freq=1;
p.Amp=1;
p.noise=.5;

[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
R=getPermutations(p.Nresidual,p.Nmeas,p.Nperm,p.Nacro);
Y=getSimulatedData(t,p);
[X,I3,I4]=constructUtilMats(t,p);
getPermutedData(Y,R,I3,I4);
toc