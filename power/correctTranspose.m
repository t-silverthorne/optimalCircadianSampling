y1=Y(1,:)
sum(y1)
y1(I(1,:,1))
sum(y1(I(1,:,1)))

Z=Y(I);
sum(Y(I(:,:,1)),2)
%%
m=2;n=3;
A=[10 20 30; 40 50 60];
I(:,:,1)=[1 3 2; 2 3 1];
I(:,:,2)=[3 2 1; 3 1 2];
offMat=repmat((0:m-1)',1,n)*n;
I+offMat;

Ap=A';
AI=pagetranspose(Ap(pagetranspose(I+offMat)))

%%
A([3 2 1; 3+2 3+1 3+3])