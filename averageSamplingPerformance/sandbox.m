addpath('../utils')
clear
parpoole
tic
% generate parameter sets
Nparamsets=50;
nreps=50;

Nleft_list=[4 6 9 12];
Nright_list=[4 6 9 12];
nu_data=NaN(numel(Nleft_list),numel(Nright_list));
unif_data=NaN(numel(Nleft_list),numel(Nright_list));

per1=.9;
per2=.2;
parfor ii=1:numel(Nleft_list)
    nu_data_loc=NaN(1,numel(Nright_list))
    unif_data_loc=NaN(1,numel(Nright_list))
    for jj=1:numel(Nright_list)
        Nleft=Nleft_list(ii);
        Nright=Nright_list(jj);
        [nu_data_loc(jj),unif_data_loc(jj)]=get_nested_average_error(per1,per2,Nleft,Nright);
    end
    nu_data(ii,:)=nu_data_loc;
    unif_data(ii,:)=unif_data_loc;
end
save('batch_sampling_test.mat','unif_data','nu_data')
toc
function [error_nu,error_unif]=get_nested_average_error(per1,per2,Nleft,Nright)
Nparamsets=1;
nreps=1;
circParams=rand(Nparamsets,4);
circParams(:,1:2)= 10.^(2*circParams(:,1:2)-1); % amplitudes 
circParams(:,3:4)=2*pi*circParams(:,3:4);   % acrophases 
linParams=get_linear_params(circParams);    % convert using trig identity


% measurement times are on the interval [0,1]
Ntimes=Nleft+Nright;

mt1=linspace(0,1/3,Nleft+1);
mt1=mt1(1:end-1);
mt2=linspace(1/3,1,Nright+1);
mt2=mt2(1:end-1);
mt_nu=[mt1 mt2]; % construct non-uniform grid 
mt_unif=linspace(0,1,Ntimes+1); % construct uniform grid 
mt_unif=mt_unif(1:end-1);

% generate data

a1=linParams(:,1);
a2=linParams(:,2);
a3=linParams(:,3);
a4=linParams(:,4);
get_Xdat = @(zts) a1.*sin(2*pi*zts/per1)+a2.*cos(2*pi*zts/per1) + a3.*sin(2*pi*zts/per2)+a4.*cos(2*pi*zts/per2) + randn(Nparamsets,numel(zts),nreps);
Xdat_unif=get_Xdat(mt_unif);
Xdat_nu=get_Xdat(mt_nu);

Xdat_unif=permute(Xdat_unif,[3 2 1]);
Xdat_nu=permute(Xdat_nu,[3 2 1]);

error_nu  =NaN(1,Nparamsets);
error_unif=NaN(1,Nparamsets);

for ii=1:Nparamsets
    error_nu(ii)=get_mean_period_error(mt_nu,Xdat_nu(:,:,ii),per1,per2);
    error_unif(ii)=get_mean_period_error(mt_unif,Xdat_unif(:,:,ii),per1,per2);
end
error_nu=mean(error_nu);
error_unif=mean(error_unif);
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

function mean_period_error=get_mean_period_error(zts,Y,per1,per2)


dp=.1;
pergrid=.1:dp:1;
for ii=pergrid
    for jj=ii+dp:dp:1
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
gof=cell(1,size(Y,1));
perSmall_vals=NaN(1,size(Y,1));
perBig_vals=NaN(1,size(Y,1));

tSmall=min(per1,per2);
tBig=max(per1,per2);

for ii=1:size(Y,1)
    [fit_result,gof{ii}]=nonlinfit(zts,Y(ii,:) ,ii_best,jj_best);
    % make per1 the smaller one
    per1loc=fit_result.per1; 
    per2loc=fit_result.per2;
    perSmall_vals(ii)=min(per1loc,per2loc);
    perBig_vals(ii)=max(per1loc,per2loc);
end

mean_period_error= mean(0.5*abs((perSmall_vals-tSmall)/tSmall) + abs((perBig_vals-tBig)/tBig));
end