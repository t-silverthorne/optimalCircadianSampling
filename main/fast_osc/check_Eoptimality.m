addpath('../utils/')

%%
N=8;
%t=linspace(0,1,N+1);
t=sort(rand(1,N));
t=t(1:end-1);
p.freq=2;
p.Amp=1;
X=constructXreduced(t,p);

fprintf('\nMin sv  of design mat:   %f\n',min(svd(X))^2)
fprintf('Min sv of Fisher mat:    %f\n',min(eig(X'*X)))
%fprintf('Min eig of Fisher mat:   %f\n',sqrt(min(eig(X'*X))))

[~,Lmin]=fminbnd(@(x) get_lambda(t,p,x),0,2*pi);
fprintf('Non centrality:          %f\n',Lmin)

function lambdapap=get_lambda(t,p,acro)
csq       = cos(2*pi*p.freq*t-acro);
lambdapap = p.Amp^2*csq*csq';
end

function X = constructXreduced(t,p)
% build design matrix X
x1=sin(2*pi*p.freq*t);
x2=cos(2*pi*p.freq*t);
X= [x1' x2'];
end
