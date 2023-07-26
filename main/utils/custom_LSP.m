function Pomega=custom_LSP(t,y,omegas)
% get least squares periodogram
N=length(t);
t=reshape(t,1,1,N);   
y=reshape(y,1,1,N);

omegas=reshape(omegas,1,1,1,length(omegas));
R=sum([cos(omegas.*t);sin(omegas.*t)].*[cos(omegas.*t) sin(omegas.*t)],3);
r=sum([cos(omegas.*t);sin(omegas.*t)].*y,3);
Pomega=(1/N)*reshape(pagemtimes(pagetranspose(r),pagemldivide(R,r)),length(omegas),1);
end
