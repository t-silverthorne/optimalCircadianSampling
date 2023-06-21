% fast implementation of Fisher-Yates shuffle
d=7;
n=10*10*10;
%Y=randn(n,d);
Y=repmat(1:d,n,1);
%Y=repmat(10:10:d*10,n,1); % construct test matrix
%Y=Y+repmat((1:n)',1,d);


% permutation
tic
for ii=1:d-1
    avec=ii-1+randi(d+1-ii,[1 n]);
    Ytemp=Y(sub2ind(size(Y),1:n,repmat(ii,1,n)));
    Y(sub2ind(size(Y),1:n,repmat(ii,1,n)))=Y(sub2ind(size(Y),1:n,avec));
    Y(sub2ind(size(Y),1:n,avec))=Ytemp;
end
toc

sum(pagetranspose(reshape(Y',d,10,10,10)),2)