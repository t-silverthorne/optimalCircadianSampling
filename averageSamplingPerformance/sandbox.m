addpath('../utils')
% convert acrophase and amplitude to linear parameters
Nparamsets=10;
nreps=5;
circParams=rand(Nparamsets,4);
circParams(:,1)= 10.^(2*circParams(:,1)-1); % amplitudes [10^1, 10^3]
circParams(:,3:4)=2*pi*circParams(:,3:4); % acrophases on [0,2*pi]
linParams=get_linear_params(circParams);


% measurement times are on the interval [0,1]
Nleft=8;
Nright=10;
Ntimes=Nleft+Nright;

mt1=linspace(0,1/3,Nleft+1);
mt1=mt1(1:end-1);
mt2=linspace(1/3,1,Nright+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2];

% construct uniform grid 
mt_unif=linspace(0,1,Ntimes+1);
mt_unif=mt_unif(1:end-1);

per1=1;
per2=.1;
a1=linParams(:,1);
a2=linParams(:,2);
a3=linParams(:,3);
a4=linParams(:,4);
regularize=true;

get_Xdat = @(zts) a1.*sin(2*pi*zts/per1)+a2.*cos(2*pi*zts/per1) + a3.*sin(2*pi*zts/per2)+a4.*cos(2*pi*zts/per2) + .1*randn(Nparamsets,numel(zts),nreps);
Xdat_unif=get_Xdat(mt_unif);
Xdat_nu=get_Xdat(mt_nu);

Xdat_unif=permute(Xdat_unif,[3 2 1]);
Xdat_nu=permute(Xdat_nu,[3 2 1]);
plot(mt_unif,Xdat_unif(1,:,1))


% testing that linear part works
zts=mt_nu;
Y=Xdat_nu(:,:,1);
ii=per1;
jj=per2
x1=sin(2*pi*zts/ii);
x2=cos(2*pi*zts/ii);
x0=ones(1,numel(zts));
x3=sin(2*pi*zts/jj);
x4=cos(2*pi*zts/jj);
X= [x0' x1' x2' x3' x4'];
betas=(X'*X)\X'*Y'
mu=mean(betas,2)'
mu(2:end)
linParams(1,:)
fits=(X*betas)'
SSres=sum((fits-Y).^2,2)
%%


dp=.1;
pergrid=.1:dp:1;

for ii=pergrid
    for jj=ii+dp:1
        x1=sin(2*pi*zts/ii);
        x2=cos(2*pi*zts/ii);
        x0=ones(1,numel(zts));
        x3=sin(2*pi*zts/jj);
        x4=cos(2*pi*zts/jj);
        X= [x0' x1' x2' x3' x4'];
        
        % do linear regression 
        betas=(X'*X)\X'*Y'; 
        fits=(X*betas)';
        SSres=sum((fits-Y).^2,2);
        SSresavg=mean(SSres);
        
        if ii==min(pergrid) && jj==min(pergrid)+dp
            ii_best=ii;
            jj_best=jj;
            best_SSresavg=SSresavg;
        elseif SSresavg < best_SSresavg 
            ii_best=ii;
            jj_best=jj;
            best_SSresavg=SSresavg;
        end
    end
end

ii_best
jj_best
nonlinfit(zts,Y(1,:) ,ii_best,jj_best)
linParams(1,:)
%%

res_nu=cell(nreps);
res_unif=cell(nreps);




for ii=1:nreps
    [res_nu{ii},gof_nu{ii}]=nonlinfit_grid_fast(mt_nu, Xdat_nu(:,:,ii),regularize);
    [res_unif{ii},gof_unif{ii}]=nonlinfit_grid_fast(mt_unif, Xdat_unif(:,:,ii),regularize);
end

function pout=get_linear_params(pin) 
% a1 = sin coeff period 1
% a2 = cos coeff period 1
% a1 = sin coeff period 2
% a2 = cos coeff period 2
A1=pin(:,1);
A2=pin(:,2);
acro1=pin(:,3);
acro2=pin(:,4);

a1= A1.*cos(acro1);
a2=-A1.*sin(acro1);
a3= A2.*cos(acro2);
a4=-A2.*sin(acro2);

pout=[a1 a2 a3 a4];
end