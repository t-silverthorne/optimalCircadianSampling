dd = @(i,j) (i==j);% could be sped up using kroneckerDelta from v2023

%% want to evaluate E_pi[ y' pi' A pi y ]
N=4;
A=rand(N,N);
A=eye(N);%A*A';
y=(A(:,1)+A(:,2));%rand(N,1);
rho1=factorial(N)/factorial(N-1);
rho2=factorial(N)/factorial(N-2);
rho3=factorial(N)/factorial(N-2);
rho4=factorial(N)/factorial(N-3);
rho5=factorial(N)/factorial(N-4);

% only one fiber
t1=@(i,j,k,l,m,n,r,s) dd(i,k)*dd(k,m)*dd(m,r)*dd(j,l)*dd(l,n)*dd(n,s)/rho1;

% fiber of 3 and fiber of 1
t2=@(i,j,k,l,m,n,r,s) ( (dd(i,k)*dd(k,m)+dd(i,k)*dd(k,r))*(1-dd(m,r))+(dd(i,m)*dd(m,r)+dd(k,m)*dd(m,r))*(1-dd(i,k)) )* ... 
                              ((dd(j,l)*dd(l,n)+dd(j,l)*dd(l,s))*(1-dd(n,s))+(dd(j,n)*dd(n,s)+dd(l,n)*dd(n,s))*(1-dd(j,l)))/rho2;

% two fibers of 2
t3=@(i,j,k,l,m,n,r,s) ((dd(i,k)*dd(m,r)+dd(i,r)*dd(k,m))*(1-dd(i,m))+dd(i,m)*dd(k,r)*(1-dd(i,k)))* ...
                              ((dd(j,l)*dd(n,s)+dd(j,s)*dd(l,n))*(1-dd(j,n))+dd(j,n)*dd(l,s)*(1-dd(j,l)))/rho3;

% fiber of two and two fibers of 1
t4=@(i,j,k,l,m,n,r,s) ((dd(i,k)+dd(k,m)-2*dd(i,k)*dd(k,m))*(1-dd(m,r))*(1-dd(r,i)) + ...
                          (dd(m,r)+dd(r,i)-2*dd(m,r)*dd(r,i))*(1-dd(i,k))*(1-dd(k,m)) + ...
                          (dd(i,m)+dd(r,k)-2*dd(i,m)*dd(r,k))*(1-dd(i,r))*(1-dd(m,k)) )* ...
                            ((dd(j,l)+dd(l,n)-2*dd(j,l)*dd(l,n))*(1-dd(n,s))*(1-dd(s,j)) + ...
                                (dd(n,s)+dd(s,j)-2*dd(n,s)*dd(s,j))*(1-dd(j,l))*(1-dd(l,n)) + ...
                                (dd(j,n)+dd(s,l)-2*dd(j,n)*dd(s,l))*(1-dd(j,s))*(1-dd(n,l)) )/rho4;

% four distinct fibers
t5=@(i,j,k,l,m,n,r,s) (1-dd(i,k))*(1-dd(i,m))*(1-dd(i,r))*(1-dd(k,m))*(1-dd(k,r))*(1-dd(m,r))* ...
                        (1-dd(j,l))*(1-dd(j,n))*(1-dd(j,s))*(1-dd(l,n))*(1-dd(l,s))*(1-dd(n,s))/rho5;
resExact=0;
t1c=0;t2c=0;t3c=0;t4c=0;t5c=0;
countoff=0;
for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                for m=1:N
                    for n=1:N
                        for r=1:N
                            for s=1:N
                                if t1(i,j,k,l,m,n,r,s)>0
                                    t1c=t1c+1;
                                end
                                if t2(i,j,k,l,m,n,r,s)>0
                                    t2c=t2c+1;
                                end
                                if t3(i,j,k,l,m,n,r,s)>0
                                    t3c=t3c+1;
                                end
                                if t4(i,j,k,l,m,n,r,s)>0
                                    t4c=t4c+1;
                                end
                                if t5(i,j,k,l,m,n,r,s)>0
                                    t5c=t5c+1;
                                end


                                term=y(j)*A(i,k)*y(l)*y(n)*A(m,r)*y(s)*( t1(i,j,k,l,m,n,r,s) + ...
                                        t2(i,j,k,l,m,n,r,s) + t3(i,j,k,l,m,n,r,s) + t4(i,j,k,l,m,n,r,s) + ...
                                        t5(i,j,k,l,m,n,r,s) );
                                if (j~=l | n~=s) & term>0 & t3(i,j,k,l,m,n,r,s)>0
                                    fprintf('%d %d %d %d %d %d %d %d \n',i,k,m,r,j,l,n,s)
                                    countoff=countoff+term;
                                end
                                if rho1*t1(i,j,k,l,m,n,r,s)+rho2*t2(i,j,k,l,m,n,r,s)+rho3*t3(i,j,k,l,m,n,r,s)+...
                                     rho4*t4(i,j,k,l,m,n,r,s)+rho5*t5(i,j,k,l,m,n,r,s) > 1
                                     fprintf('%d %d %d %d %d %d %d %d \n',i,k,m,r,j,l,n,s)
                                end
%                                 if abs(term)>0
%                                     bool=(i==k)*(m==r)*(j==1 | j==2)*(l==1 | l==2)*(n==1 | n==2)*(s==1 | s==2);
%                                     fprintf('%d %d %d %d %d %d %d %d        %d \n',i,k,m,r,j,l,n,s,bool)
%                                 end
                                resExact=resExact+term;
                            end
                        end
                    end
                end
            end
        end
    end
end

% t1c
% t2c==(nchoosek(4,2)*nchoosek(4,3)*2)^2
% t3c/36/36
% t4c==(nchoosek(4,3)*nchoosek(4,2)*factorial(3))^2
% t5c/factorial(4)^2
% 

countoff
resExact


resNum  =0;
P=perms(1:N);
pimat=eye(N);
for ii=1:size(P,1)
    piloc=pimat(:,P(ii,:));
    resNum=resNum + (y'*piloc'*A*piloc*y)^2;
end
resNum=resNum/factorial(N);

resNum

