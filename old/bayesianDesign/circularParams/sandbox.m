[unif,~,nu,~]=estimateBayesianFIMdetCirc(2,6,0,1/3,'sdev');
mean(unif)
mean(nu)
%% setup
update_rule='MCMC' % options: MCMC or regression
NL=3;
NR=8;
tauA=0;
tauB=1/3;
[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

param.deltaT=.1;
param.A1max=100;
param.A2max=100;
ptc=num2cell(paramTrue);


%% simulate experiment
biCosinor=@(t,A1,phi1,T1,A2,phi2,T2) A1*cos(2*pi*t/T1-phi1)+A2*cos(2*pi*t/T2-phi2)+randn(1,numel(t));
paramTrue=[1,pi,.5,1,0,.8];
ptc=num2cell(paramTrue); % true parameters stored as cell

Yobs_unif=biCosinor(mt_unif,ptc{:});
Yobs_nu=biCosinor(mt_nu,ptc{:});



%% update prior or do regression
switch update_rule
    case 'MCMC'
        posterior_unif=runMCMC(10,@(A1,phi1,T1,A2,phi2,T2,Yobs,tobs,param) ...
            almostUnifImproperPrior(A1,phi1,T1,A2,phi2,T2,param)*...
            biharmonicLikelihood(A1,phi1,T1,A2,phi2,T2,Yobs,tobs) ...
            ,Yobs_unif,mt_unif,param);
        posterior_nu=runMCMC(10,@(A1,phi1,T1,A2,phi2,T2,Yobs,tobs,param) ...
            almostUnifImproperPrior(A1,phi1,T1,A2,phi2,T2,param)*...
            biharmonicLikelihood(A1,phi1,T1,A2,phi2,T2,Yobs,tobs) ...
            ,Yobs_unif,mt_unif,param);
end

%% optimize non-uniform strategy

function p=almostUnifImproperPrior(A1,phi1,T1,A2,phi2,T2,param)
% return the modified uniform prior up to normalization constant
if A1>param.A1max || A1<0 || A2>param.A2max  || A2<0 ||  ...
    phi1<0 || phi1>2*pi || phi2<0 || phi2>2*pi ||  ...
    T1<0 || T1>1 || T2<0 || T2>1|| ...
    abs(T1-T2)<param.deltaT
    
    p=0;
else
    p=1;
end
end

function p=biharmonicLikelihood(A1,phi1,T1,A2,phi2,T2,Yobs,tobs)
p=prod(normpdf(A1*cos(2*pi*tobs/T1-phi1)+A2*cos(2*pi*tobs/T2-phi2)-Yobs,0,1));
end

function Sampvec=runMCMC(Nsamp,p,Yobs,tobs,param)
T=1e2; % length of Markov chain
Sampvec=NaN(Nsamp,6);
for N=1:Nsamp
    Xt=rand(1,6);
    for t=1:T
        Y=randn(1,6)+Xt;
        ntcY=num2cell(Y); ntcXt=num2cell(Xt);
        alpha=min( p(ntcY{:},Yobs,tobs,param)*normpdfvec(Xt,Y)/ p(ntcXt{:},Yobs,tobs,param)/normpdfvec(Y,Xt),1 );
        if rand<alpha
            Xt=Y;
        end
    end
    Sampvec(N,:)=Xt;
end
end

function p = normpdfvec(X,Y)
p=prod(arrayfun(@(ind) normpdf(X(ind),Y(ind)),1:numel(X)));
end


% function Sampvec=almostUnifPrior(Nsamp,param)
% Sampvec=NaN(Nsamp,6); % parameter matrix
% Sampvec(:,[1 4])=10*rand(Nsamp,2); % get amplitudes
% Sampvec(:,[2 5])= 2*pi*rand(Nsamp,2); % get acrophases
% accepted=NaN(Nsamp,2); 
% % final two need rejection sampling
% ind=1;
% while ind<Nsamp
%     samp=rand(Nsamp,2);
%     acceptedloc=samp(abs(samp(:,1)-samp(:,2))>param.deltaT,:);
%     accepted(ind:ind+size(acceptedloc,1)-1,:)=acceptedloc;
%     ind=ind+size(acceptedloc,1)-1;
% end
% Sampvec(:,[3 6])=accepted(1:Nsamp,:); % get periods
% end

