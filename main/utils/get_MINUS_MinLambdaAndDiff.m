function [lam,DlamDt] = get_MINUS_MinLambdaAndDiff(mt,freq)
% return min lambda value and its derivative wrt measurement times
lam=-1*getMinLambda(mt,freq);
if nargout >1
  DlamDt=-1*diffMinLambdaDt(mt,freq);
end
end
