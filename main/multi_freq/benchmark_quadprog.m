close all
clear
addpath('../utils/')

f0    = 3.8;        % frequency in data
Amp   = 1.5;     % amp
Nmeas = 8;

use_eq_constraints=true;
% measurement times (coarse grid)
t     = linspace(0,1,Nmeas+1); 
t     = t(1:end-1);
t=t';

% frequency grid
omegas= 2*pi*(0:.1:7);  
omegas= omegas(2:end-1); % useful to avoid end frequencies

% time grid for optimization (fine grid)
tau=(0:.005:1)'; % candidate measurement times
tau=tau(1:end-1);
mu0= ones(length(t),1); % existing times not being varied
mu = ones(length(tau),1)/length(tau); % weights for candidate meas times
t3d=reshape(t,[1,1,length(t)]);

% data
y=Amp*cos(2*pi*f0*t-0.2)+randn(length(t),1);

% estimate LSP
Pomegas=custom_LSP(t,y,omegas);




%% construct quadprog
% constant term (useful for comparison but redundant for optimization)
a=0;
for ii=1:length(omegas)
    for jj=1:length(omegas)
        a=a+sum(Pomegas(ii)*Pomegas(jj)*exp(1j*(omegas(ii)-omegas(jj))*(t-t')),'all');
    end
end
a=real(a);

% linear term
f=arrayfun(@(ind) sum((Pomegas*Pomegas').*exp(1j*(omegas-omegas').*(tau(ind)-t3d)),'all'),(1:length(mu))');
f=f+conj(f);

% quadratic term
Htau=zeros(length(tau),length(tau));
for ii=1:length(omegas)
    for jj=1:length(omegas)
        Htau=Htau+Pomegas(ii)*Pomegas(jj)*exp(1j*(omegas(ii)-omegas(jj))*(tau-tau'));
    end
end
Htau=real(Htau);

% run quadprog


% without equality constraints

if use_eq_constraints
    Aeq=ones(1,length(mu)); % require normalization of distribution
    beq=1;
    [muopt,copt]=quadprog(2*Htau,f,[],[],...
          Aeq,beq, ...
          zeros(length(mu),1),ones(length(mu),1),mu);
else
    [muopt,copt]=quadprog(2*Htau,f,[],[],...
          [],[], ...
          zeros(length(mu),1),ones(length(mu),1),mu);
end

% get candidate measurement times
[~,opt_inds]=maxk(muopt,Nmeas);

topt=transpose(tau(opt_inds));

% plot result
tiledlayout(1,2)
nexttile(1)
xline(topt)
nexttile(2)

numFgrid=2^7; % for contour plots
numAgrid=2^4;
freqvals  = linspace(1,18,numFgrid);
Ampvals   = logspace(-1,1,numAgrid);
[Agr,Fgr]=meshgrid(Ampvals,freqvals);
Pmin=arrayfun(@(Agr,Fgr) getMinPower(Agr,Fgr,topt),Agr, Fgr);
[M,c]=contourf(Agr,Fgr,Pmin,100,'LineColor','none');%,'FaceAlpha',0.
hold on % lines for indicating cross sections
ylim([freqvals(1), freqvals(end)])
set(gca,'XScale','log') % rescale x axis
xlabel('amplitude') % xlabel
set(gca,'YTickLabel',[]);


% add colorbar shared by first two plots
set(cont, 'Colormap',parula, 'CLim', [0 1])
% assign color bar to one tile 
cbtop = colorbar(cont(end));
cbtop.Label.Interpreter='latex';
set(cbtop,'TickLabelInterpreter','latex')
cbtop.Label.String='$\textrm{min}_\phi \gamma(\phi;A,f)$';
cbtop.Label.Interpreter='latex';

