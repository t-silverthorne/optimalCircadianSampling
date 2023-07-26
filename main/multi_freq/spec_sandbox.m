P=@(x) x^2;
N=1e2;
t=rand(1e2,1);
omegas=0:1:10;
P=@(x)x.^2;


t=(0:.05:1)';
f0=4
y=2.3*cos(2*pi*f0*t-0.2)+.1*randn(length(t),1);
t=reshape(t,1,1,length(t));   
y=reshape(y,1,1,length(y));
% estimate LSP
omegas=2*pi*(0:.1:7);
omegas=omegas(2:end-1)
tic
Pomega=custom_LSP(t,y,omegas);
toc
clf
tiledlayout(2,1)
nexttile(1)
plot(omegas/2/pi,Pomega,'-ok')
xline(f0)
tic
[pxx,f]=plomb(reshape(y,length(y),1),reshape(t,length(t),1), ...
                omegas/2/pi,[],'power');
toc
hold on
plot(f,pxx,'-or')
nexttile(2)
semilogy(f,abs((pxx-Pomega)),'.k')

mu=ones(length(t),1)/length(t);
tic
Jspec_wL2_slow(t,mu,omegas,Pomega);
toc
tic
Jspec_wL2(t,mu,omegas,Pomega);
toc
%% quadprog setup

tau=0:.01:1;% candidate times
t=reshape(t,[1,1,length(t)])
f=zeros(length(mu),1);

f=arrayfun(@(ind) sum((Pomega-Pomega').*exp(1j*(omegas-omegas').*(t-tau(ind))),'all'),(1:length(mu))')

%%
omega=2*pi*3.5;
R=sum([cos(omega*t);sin(omega*t)].*[cos(omega*t) sin(omega*t)],3);
r=sum([cos(omega*t);sin(omega*t)].*y);
R\r

function J = Jspec_wL2_slow(t,mu,omegas,Pomegas)
J=0;
t=reshape(t,length(t),1);
mu=reshape(mu,length(mu),1);
for p=1:length(omegas)
    for k=1:length(omegas)
        J=J+Pomegas(p)*Pomegas(k)*abs(sum(exp(1j*(omegas(p)-omegas(k))*t).*mu))^2;
    end
end
end


function J = Jspec_wL2(t,mu,omegas,Pomegas)
% weighted L2 cost function for spectral optimization
t=reshape(t,1,1,length(t));   
mu=reshape(mu,1,1,length(mu));
J1=exp(1j*(omegas-omegas').*t);
J1=sum(mu.*J1,3);
J1=J1.*conj(J1);
J2=Pomegas*Pomegas';
J=sum(J1.*J2,'all');
end

function Pomega=custom_LSP(t,y,omegas)
% get least squares periodogram
N=length(t);
omegas=reshape(omegas,1,1,1,length(omegas));
R=sum([cos(omegas.*t);sin(omegas.*t)].*[cos(omegas.*t) sin(omegas.*t)],3);
CC=arrayfun(@(ind) cond(R(:,:,ind)),1:length(omegas));
r=sum([cos(omegas.*t);sin(omegas.*t)].*y,3);
Pomega=(1/N)*reshape(pagemtimes(pagetranspose(r),pagemldivide(R,r)),length(omegas),1);
end
