function p = normpdfvec(X,Y,sdvec)
if nargin<3
    sdvec=ones(1,numel(X));
end
p=prod(arrayfun(@(ind) normpdf(X(ind),Y(ind),sdvec(ind)),1:numel(X)));
end