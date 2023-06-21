function P = getPermutations(S1,d,S3,S4)
% GETPERMUTATIONS Fast implementation of Fisher-Yates algorithm for
% generating permutations on symbols {1,...,d}. This generates S1xS3xS4
% permutations of d symbols.
% At the end of the implementation, the array of permutations is reshaped
% to have size S1 x d x S3 x S4, so that it is compatible with the
% getPermutedData function.
cheap_perm=false;
if cheap_perm
    S4old=S4;
    S4=1;
end
n=S1*S3*S4;
P=repmat(1:d,n,1);

P=P';
for ii=1:d-1
    avec=ii-1+randi(d+1-ii,[n 1]);
    for jj=1:n
        Ytemp=P(ii,jj);
        P(ii,jj)=P(avec(jj),jj);
        P(avec(jj),jj)=Ytemp;
    end
end
P=P';
P=pagetranspose(reshape(P',d,S1,S3,S4));
if cheap_perm
    P=repmat(P,1,1,1,S4old);
end
end
