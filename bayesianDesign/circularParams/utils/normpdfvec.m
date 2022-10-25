function p = normpdfvec(X,Y)
p=prod(arrayfun(@(ind) normpdf(X(ind),Y(ind)),1:numel(X)));
end