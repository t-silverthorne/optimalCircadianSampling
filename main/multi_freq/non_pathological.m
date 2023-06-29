% checking claim that R -> N/2 I for large N
addpath('../utils/')
tiledlayout(2,1)
nreps=1e1;
Nmax=1e2;
%% uniform random
clf
for N=1:Nmax
    for k=1:nreps
        t=sort(rand(1,N));
        f=2;
        omega=2*pi*f;
        t=t/f;
        CC=cos(omega*t)*cos(omega*t)';
        CS=cos(omega*t)*sin(omega*t)';
        SS=sin(omega*t)*sin(omega*t)';
        
        nexttile(1)
        plot(N,norm([CC CS; CS SS] -N/2*eye(2)),'.k')
        hold on
        nexttile(2)
        semilogy(N,cond([CC CS; CS SS]),'.k')
        ylim([1,1e3])
        hold on
   end
end

%% 2-uniform
clf
Nmax=1e2
for N=1:Nmax
        [~,t]=getSamplingSchedules(N-floor(N/3),floor(N/3),0,0.3)
        f=1;
        omega=2*pi*f;
        t=t/f;
        CC=cos(omega*t)*cos(omega*t)';
        CS=cos(omega*t)*sin(omega*t)';
        SS=sin(omega*t)*sin(omega*t)';
        
        nexttile(1)
        plot(N,norm([CC CS; CS SS] -N/2*eye(2))/N,'.k')
        hold on
        nexttile(2)
        semilogy(N,cond([CC CS; CS SS]),'.k')
        hold on
end