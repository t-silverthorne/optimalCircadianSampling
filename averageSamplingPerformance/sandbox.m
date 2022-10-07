addpath('../utils')
clear
tic
% generate parameter sets

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
save('batch_sampling_test_500_250.mat','unif_data','nu_data')
toc
function [error_nu,error_unif]=get_nested_average_error(per1,per2,Nleft,Nright)
Nparamsets=500;
nreps=250;
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

sdev_nu=Inf;
sdev_unif=Inf;
while sdev_nu || sdev_unif > 1e-1
    error_nu(ii)=get_mean_period_error(mt_nu,Xdat_nu(:,:,ii),per1,per2);
    error_unif(ii)=get_mean_period_error(mt_unif,Xdat_unif(:,:,ii),per1,per2);
		sdev_nu=std(error_nu);
		sdev_unif=std(error_unif);
end
error_nu=mean(error_nu);
error_unif=mean(error_unif);
end

