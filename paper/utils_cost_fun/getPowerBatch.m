function [pwr_mat,amp_mat,acro_mat]=getPowerBatch(t,p,I3,I4)
%%%%%%%%%%%%%%
%GETPOWERBATCH wrapper for getPower function
% assumes that struct p has field 'Nbatch'
% returns matrices whose columns correspond to acrophases and rows
% correspond to batches
%%%%%%%%%%%%%%
pwr_mat=[];amp_mat=[];acro_mat=[];
for ii=1:p.Nbatch
    [pwr_est,amp_est,acro_est]=getPower(t,p,I3,I4);
    pwr_est  = reshape(pwr_est,1,p.Nacro);
    amp_est  = reshape(amp_est,[p.Nresidual,p.Nacro]);
    acro_est = reshape(acro_est,[p.Nresidual,p.Nacro]);
    if isempty(acro_mat)
        pwr_mat  = pwr_est;
        acro_mat = acro_est;
        amp_mat  = amp_est;
    else
        pwr_mat  = [pwr_mat;pwr_est];
        acro_mat = [acro_mat;acro_est];
        amp_mat  = [amp_mat;amp_est];
    end
end
end