% run nonlinear regression with different initial guesses

clf
clear
addpath('../utils/')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

warning('off','MATLAB:nearlySingularMatrix')
pctRunOnAll warning off

method=2; % method that gave performance_difference
Nleft=8;
Nright=2;
regularize=true; % only consider periods greater than 2hrs (otherwise aliasing occurs)
for outer_index=1:4
	disp(outer_index)
    switch method
        case 1
            Ntimes=14; 
            frq1=1;  % samples per hour in first region
            frq2=.3; % sph in second region
            f= [frq1 frq2 -1];
            A=-[frq1 frq2  0; zeros(2,3)];
            b=-[Ntimes-2; zeros(2,1)];
            Aeq=[1 1 0; 0 0 1; 0 0 0];
            beq=[24; Ntimes;0];
            lb=[1 1 1];
            nui=intlinprog(f,1:3,A,b,Aeq,beq,lb);
            Ntimes=ceil(nui(1:2)'*[frq1 frq2 ]');
            
            mt1=linspace(0,nui(1),frq1*nui(1)+1);
            mt1=mt1(1:end-1);
            mt2=linspace(nui(1),nui(1)+nui(2),ceil(frq2*nui(2))+1);
            mt2=mt2(1:end-1);
            mt_nu=[mt1 mt2];
            
            % construct uniform grid 
            mt_unif=linspace(0,24,Ntimes+1);
            mt_unif=mt_unif(1:end-1);
            
            % sampling options
            Nsamples=1; % number of samples at each time point
            nreps=1e2;  % number of times to repeat experiment
        % better setup for looping over indices
        case {2,3} 
            Nleft=Nleft+2;
            Nright=Nright+2;
            Ntimes=Nleft+Nright;
            
            mt1=linspace(0,8,Nleft+1);
            mt1=mt1(1:end-1);
            mt2=linspace(8,24,Nright+1);
            mt2=mt2(1:end-1);
            mt_nu=[mt1 mt2];
            
            % construct uniform grid 
            mt_unif=linspace(0,24,Ntimes+1);
            mt_unif=mt_unif(1:end-1);
            
            % sampling options
            Nsamples=1; % number of samples at each time point
            nreps=1e3;  % number of times to repeat experiment
    end
    per1=12;
    per2=4;
    zts_unif = repmat(mt_unif,1,Nsamples);
    zts_nu   = repmat(mt_nu,1,Nsamples);
    a0=0; a1=0; a2=.5; a3=.25; a4=.25;
    sig=0.1; % noise level
    
    avec=[a0 a1 a2 a3 a4 per1 per2]; % for plotting
    
    % for generating data
    get_Xdat = @(zts) a0+a1*sin(2*pi*zts/per1)+a2*cos(2*pi*zts/per1) + a3*sin(2*pi*zts/per2)+a4*cos(2*pi*zts/per2) + sig*randn(nreps,numel(zts));
    
    Xdat_unif=get_Xdat(zts_unif); % sample on uniform grid
    Xdat_nu=get_Xdat(zts_nu); % sample on non-uniform grid
    

     
%     param.beta0=a0;
%     param.beta1=a1;
%     param.beta2=a2;
%     param.beta3=a3;
%     param.beta4=a4;
%     param.T1=per1/24;
%     param.T2=per2/24;
%     log(det(get_FIM_biharmonic_nonlin(mt_nu',ones(numel(mt_nu),1)/numel(mt_nu),param)))
     


    res_nu=cell(nreps,1);
    res_unif=cell(nreps,1);
    gof_nu=cell(nreps,1);
    gof_unif=cell(nreps,1);
    
%     clf
%     scatter(zts_unif,Xdat_unif(1,:))
%     hold on
%     scatter(zts_nu,Xdat_nu(1,:))
%     zts=0:.01:24;
%     plot(zts,a0+a1*sin(2*pi*zts/per1)+a2*cos(2*pi*zts/per1) + a3*sin(2*pi*zts/per2)+a4*cos(2*pi*zts/per2))
%     hold off
    
    Xdat_unif=parallel.pool.Constant(Xdat_unif);
    Xdat_nu=parallel.pool.Constant(Xdat_nu);
    
    
    switch method
        case {1,2}
            tic
            parfor ii=1:nreps
                [res_nu{ii},gof_nu{ii}]=nonlinfit_grid_fast(zts_nu, Xdat_nu.Value(ii,:),regularize);
                [res_unif{ii},gof_unif{ii}]=nonlinfit_grid_fast(zts_unif, Xdat_unif.Value(ii,:),regularize);
            end
            toc
        case 3
            % run first sample and generate an initial guess at period
            [res_nu{1},gof_nu{1},per1_init_nu,per2_init_nu]=nonlinfit_grid(zts_nu, Xdat_nu.Value(1,:));
            [res_unif{1},gof_unif{1},per1_init_unif,per2_init_unif]=nonlinfit_grid(zts_unif, Xdat_unif.Value(1,:));

            % reuse initial period guess
            parfor ii=2:nreps
                [res_nu{ii},gof_nu{ii}]=nonlinfit(zts_nu, Xdat_nu.Value(ii,:),per1_init_nu,per2_init_nu);
                [res_unif{ii},gof_unif{ii}]=nonlinfit(zts_unif, Xdat_unif.Value(ii,:),per1_init_unif,per2_init_unif);
            end
            
    end
    switch method
        case 1
            save(strcat( ...
            'method_', ...
            num2str(method),...
            '_Nsamp_',...
            num2str(Nsamples),... 
            '_nreps_',...
            num2str(nreps),...
            '_per1_',...
            num2str(per1),...
            '_per2_',...
            num2str(per2), ...
            '_Ntimes_',...
            num2str(Ntimes)))
        case {2,3}
            save(strcat( ...
            'final_method', ...
            num2str(method),...
            '_regularize_',...
            num2str(regularize),...
            '_Nsamp_',...
            num2str(Nsamples),... 
            '_nreps_',...
            num2str(nreps),...
            '_per1_',...
            num2str(per1),...
            '_per2_',...
            num2str(per2), ...
            '_Nleft_',...
            num2str(Nleft), ...
            '_Nright_',...
            num2str(Nright)))
    end
end
warning('on','MATLAB:nearlySingularMatrix')
pctRunOnAll warning on
