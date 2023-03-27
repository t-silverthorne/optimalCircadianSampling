function P = fisher_yates_rank4(S1,d,S3,S4)
% DESCRIPTION:
% Fast implementation of Fisher-Yates algorithm for generating permutations
% on symbols {1,...,d}.
%
% At the end of the implementation, the matrix of permutations is reshaped
% to have size S1 x d x S3 x S4, so that it is compatible with the
% get_permuted_Y_rank4 function.
n=S1*S3*S4;
P=repmat(1:d,n,1);

% permutation

for ii=1:d-1
    avec=ii-1+randi(d+1-ii,[1 n]);
    Ytemp=P(sub2ind(size(P),1:n,repmat(ii,1,n)));
    P(sub2ind(size(P),1:n,repmat(ii,1,n)))=P(sub2ind(size(P),1:n,avec));
    P(sub2ind(size(P),1:n,avec))=Ytemp;
end

P=pagetranspose(reshape(P',d,S1,S3,S4));
end

