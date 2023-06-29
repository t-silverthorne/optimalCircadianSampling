dd = @(i,j) (i==j);% could be sped up using kroneckerDelta from v2023

% want to evaluate E_pi[ y' pi' A pi y ]
N=6;
A=rand(N,N);
%E=eye(N)
y=rand(N,1);
rho1=factorial(N)/factorial(N-1);
rho2=factorial(N)/factorial(N-2);
rho3=factorial(N)/factorial(N-2);
rho4=factorial(N)/factorial(N-3);
rho5=factorial(N)/factorial(N-4);

% only one fiber
t1=@(i,j,k,l,m,n,r,s) dd(i,k)*dd(k,m)*dd(m,r)*dd(j,l)*dd(l,n)*dd(n,s)/rho1;

% fiber of 3 and fiber of 1
t2=@(i,j,k,l,m,n,r,s) ( ( dd(i,k)*dd(k,m)*(1-dd(m,r)) )*( dd(j,l)*dd(l,n)*(1-dd(n,s)) ) + ...
                        ( dd(i,k)*dd(k,r)*(1-dd(m,r)) )*( dd(j,l)*dd(l,s)*(1-dd(n,s)) ) + ...
                        ( dd(i,m)*dd(m,r)*(1-dd(i,k)) )*( dd(j,n)*dd(n,s)*(1-dd(j,l)) ) + ...
                        ( dd(k,m)*dd(m,r)*(1-dd(i,k)) )*( dd(l,n)*dd(n,s)*(1-dd(j,l)) ))/rho2;

% two fibers of 2
t3=@(i,j,k,l,m,n,r,s) ( ( dd(i,k)*dd(m,r)*(1-dd(i,m)) )*( dd(j,l)*dd(n,s)*(1-dd(j,n)) ) + ...
                        ( dd(i,m)*dd(k,r)*(1-dd(i,r)) )*( dd(j,n)*dd(l,s)*(1-dd(j,s)) ) + ...
                        ( dd(i,r)*dd(k,m)*(1-dd(i,k)) )*( dd(j,s)*dd(l,n)*(1-dd(j,l)) ))/rho3;

% fiber of two and two fibers of 1
t4=@(i,j,k,l,m,n,r,s) (( dd(i,k)*(1-dd(k,m))*(1-dd(m,r))*(1-dd(r,i)) )*( dd(j,l)*(1-dd(l,n))*(1-dd(n,s))*(1-dd(s,j)) ) + ... 
                       ( dd(k,m)*(1-dd(m,r))*(1-dd(r,i))*(1-dd(i,k)) )*( dd(l,n)*(1-dd(n,s))*(1-dd(s,j))*(1-dd(j,l)) ) + ... 
                       ( dd(m,r)*(1-dd(r,i))*(1-dd(i,k))*(1-dd(k,m)) )*( dd(n,s)*(1-dd(s,j))*(1-dd(j,l))*(1-dd(l,n)) ) + ... 
                       ( dd(r,i)*(1-dd(i,k))*(1-dd(k,m))*(1-dd(m,r)) )*( dd(s,j)*(1-dd(j,l))*(1-dd(l,n))*(1-dd(n,s)) ) + ... 
                       ( dd(i,m)*(1-dd(i,r))*(1-dd(i,k))*(1-dd(r,k)) )*( dd(j,n)*(1-dd(j,s))*(1-dd(j,l))*(1-dd(s,l)) ) + ... 
                       ( dd(r,k)*(1-dd(i,m))*(1-dd(i,r))*(1-dd(m,k)) )*( dd(s,l)*(1-dd(j,n))*(1-dd(j,s))*(1-dd(n,l)) ))/rho4;

% four distinct fibers
t5=@(i,j,k,l,m,n,r,s) ((1-dd(i,k))*(1-dd(i,m))*(1-dd(i,r))*(1-dd(k,m))*(1-dd(k,r))*(1-dd(m,r))* ...
                        (1-dd(j,l))*(1-dd(j,n))*(1-dd(j,s))*(1-dd(l,n))*(1-dd(l,s))*(1-dd(n,s)))/rho5;

resExact=0;
t1c=0;t2c=0;t3c=0;t4c=0;t5c=0;
countoff=0;
tic
for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                for m=1:N
                    for n=1:N
                        for r=1:N
                            for s=1:N
                                term=y(j)*A(i,k)*y(l)*y(n)*A(m,r)*y(s)*( t1(i,j,k,l,m,n,r,s) + ...
                                        t2(i,j,k,l,m,n,r,s) + t3(i,j,k,l,m,n,r,s) + t4(i,j,k,l,m,n,r,s) + ...
                                        t5(i,j,k,l,m,n,r,s) );

                                resExact=resExact+term;
                            end
                        end
                    end
                end
            end
        end
    end
end
toc

tic
resExactSubInd=0;
parfor ind=1:N^8
    [i,j,k,l,m,n,r,s]=ind2sub(N*ones(1,8),ind);
    term=y(j)*A(i,k)*y(l)*y(n)*A(m,r)*y(s)*( t1(i,j,k,l,m,n,r,s) + ...
            t2(i,j,k,l,m,n,r,s) + t3(i,j,k,l,m,n,r,s) + t4(i,j,k,l,m,n,r,s) + ...
            t5(i,j,k,l,m,n,r,s) );
    resExactSubInd=resExactSubInd+term;
end
toc

resNum =0;
P=perms(1:N);
pimat=eye(N);
for ii=1:size(P,1)
    piloc=pimat(:,P(ii,:));
    resNum=resNum + (y'*piloc'*A*piloc*y)^2;
end
resNum=resNum/factorial(N);

resNum
resExact
resExactSubInd