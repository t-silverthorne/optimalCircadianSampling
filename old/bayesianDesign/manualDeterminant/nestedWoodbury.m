d=6;
nconv=7;
A=rand(d,nconv);

B=0;
for ii=1:nconv
    B=B+A(:,ii)*A(:,ii)'
end
det(A*A')

