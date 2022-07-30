clear
addpath('../utils')

Nptsvec=8:8
detvec=NaN(1,numel(Nptsvec));
for ii=1:numel(Nptsvec)
    Npts=Nptsvec(ii);
    theta_vec=linspace(0,1,Npts+1);
    theta_vec=theta_vec(1:end-1)';
    w_vec=ones(numel(theta_vec),1)/numel(theta_vec);
    detvec(ii)=-log(det(get_memory_info_matrix(theta_vec,w_vec,'tanhsq')));
end

display(detvec)
% plot(Nptsvec,detvec,'.k')
