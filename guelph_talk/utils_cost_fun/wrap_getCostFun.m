function Jvec = wrap_getCostFun(t,p,active_inds)
% WRAP_GETCOSTFUN wrapper for getCostFun
[pwr_mat,amp_mat,acro_mat]=getPowerBatch(t,p);
Jvec = getCostFun(p,pwr_mat,amp_mat,acro_mat,active_inds);
end

