function Jvec = wrap_getCostFun(t,p,I3,I4)
% WRAP_GETCOSTFUN wrapper for getCostFun
[pwr_mat,amp_mat,acro_mat]=getPowerBatch(t,p,I3,I4);
Jvec = getCostFun(p,pwr_mat,amp_mat,acro_mat);
end

