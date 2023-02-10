d=20;
N=1e5;
fprintf('\n')
tic
P=rand(d,d,N,'gpuArray');
A=rand(d,d,'gpuArray');
P=gather(P); A=gather(A);
sum(pagemtimes(pagetranspose(P),pagemtimes(A,P)),3)/N;
toc

tic
P=rand(d,d,N);
A=rand(d,d);
sum(pagemtimes(pagetranspose(P),pagemtimes(A,P)),3)/N;
toc

%% 3D array with all permutations 
d=8;
pall=perms(1:d);

Pmat=NaN(d,d,factorial(d));
I=eye(d);
for ii=1:factorial(d)
    Pmat(:,:,ii)= I(:,pall(ii,:));
end
Pmat(:,:,1)
