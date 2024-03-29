%% 
p.Nmeas     = 12;
p.Nacro     = 32;
p.Nresidual = 1e3;
p.Nperm     = 1e2;
p.noise     = 1;
p.Nbatch    = 1;
p.Amp       = 2;
p.freq      = 2.8;


Y           = rand(p.Nresidual,p.Nmeas,1,p.Nacro); % data to be permuted

% options: naive, naive_make_perms_first, naive_reuse_perms,
%permMethod       = 'naive';
permMethod        = 'FY_double_for_loop'; permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
%permMethod       = 'naive_reuse_perms'; 


doExtrapolation  = false; 

% naive
tic
switch permMethod
    case 'naive'
    % loop over all data and generate permutations one at a time (slowest method)
    Y = repmat(Y,[1 1 p.Nperm 1]);
        for kk=1:p.Nacro
            for jj=1:p.Nresidual
                for ii=1:p.Nperm
                    Y(jj,:,ii,kk)=Y(jj,randperm(p.Nmeas,p.Nmeas),ii,kk);
                end
            end
        end

    case 'FY_double_for_loop'
        Y = repmat(Y,[1 1 p.Nperm 1]);
        N1=p.Nresidual;d=p.Nmeas;N2=d;N3=p.Nperm;N4=p.Nacro;
        n=N1*N3*N4;
        I3=NaN(N1,N2,N3); % make indexing matrices
        I4=NaN(N1,N2,N3,N4);
        for ii=1:N3
            I3(:,:,ii)=ii*ones(N1,N2);
        end
        for ii=1:N4
            I4(:,:,:,ii)=ii*ones(N1,N2,N3);
        end
        
        P=repmat(1:d,n,1); % make permutation matrix using FY shuffle
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
        P=pagetranspose(reshape(P',d,N1,N3,N4));
        % reindex to act with permutation
        switch permActionMethod
            case 'index'
                Y=Y(sub2ind(size(Y),repmat((1:N1)',1,N2,N3,N4),P,repmat(I3,1,1,1,N4),I4));
            case 'matrix'
                d=N2;
                I=eye(d);
                I=repmat(I,[1 1 N3 N4]);
                ind3=repmat(ones(d,d),[1 1 N3 N4]);
                for ii=1:N3
                    ind3(:,:,ii,:)=ind3(:,:,ii,:)*ii;
                end
                ind4=repmat(ones(d,d),[1 1 N3 N4]);
                for ii=1:N4
                    ind4(:,:,:,ii)=ind4(:,:,:,ii)*ii;
                end
                Yloc=Y;
                YI=NaN(N1,N2,N3,N4); % construct in transpose form
                for ii=1:N1
                    Ploc=P(ii,:,:,:);
                    YI(ii,:,:,:)=pagemtimes(Yloc(ii,:,:,:),I(sub2ind(size(I),repmat((1:d)',[1 d N3 N4]),repmat(Ploc,[d 1 1 1]),ind3,ind4)));
                end
        end

    case 'naive_reuse_perms'
        Y = repmat(Y,[1 1 p.Nperm 1]);
        for ii=1:p.Nperm
            rp=randperm(p.Nmeas,p.Nmeas);
            I1=repmat((1:p.Nresidual)',[1 p.Nmeas,1,p.Nacro]);
            I2=repmat(rp,[p.Nresidual 1 1 p.Nacro]);
            I3=ii*ones([p.Nresidual,p.Nmeas,1,p.Nacro]);
            I4=reshape(1:p.Nacro,[1 1 1 p.Nacro]).*ones([p.Nresidual,p.Nmeas,1,p.Nacro]);
            Y(:,:,ii,:) = Y(sub2ind(size(Y),I1,I2,I3,I4));
        end
    
end
toc