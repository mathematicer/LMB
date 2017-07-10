function [ X ] = LMB(Dz, k, d_hat, mode)
%Based on the pairwise distance matrix Dz, LMB generates a d_hat-dimension 
%representation X=(x1,x2,...,xn) with x1=0. k is the neighbourhood size. l
%is the number of reduced points used to compute the orthogonal
%transformation. k>=d_hat should be satisfied for LMB to work. mode=1 for
%closed manifold, mode=2 for open manifold.
[n,~]=size(Dz);Dc=Dz;
l=k;k_LMB=k;
index(1)=1;
reduced=zeros(1,n);
reduced(1)=1;
X=zeros(d_hat,n);
front=1;rear=1;
learning=1;
while(front<=rear||front<=n)
    if(front==rear+1)
        front=front-1;
        k_LMB=k_LMB+5;
    else
        k_LMB=k;
    end
    if mode==1
        [~,IDX]=sort(Dc(index(front),:));
    else
        [~,IDX]=sort(Dz(index(front),:));
    end
    I=IDX(1:k_LMB+1);
    new=find(reduced(I)==0);
    if(isempty(new))
        front=front+1;
        continue;
    end
    J=I(new);
    Di=Dz(I,I);
    if(index(front)==1)
        [Y]=MDS(Di,d_hat);
        X(:,J)=Y(:,new)+(X(:,1)-Y(:,1))*ones(1,length(J));
        if mode==1
            Dc(index,J)=dist(X(:,index)',X(:,J));
            Dc(J,index)=dist(X(:,J)',X(:,index));
            if(length(J)~=1)
                Dc(J,J)=squareform(pdist(X(:,J)'));
            end
        end
        rear=rear+length(new);
        reduced(J)=1;
        front=front+1;
        index=[index,J];
    else
        V=knnsearch((X(:,index))',(X(:,index(front)))','K',l+1);
        close=index(V(2:l+1));
        C=X(:,index(front))*ones(1,l)-X(:,close);
        Y=MDS(Di,d_hat);
        Xr=Y(:,2:k_LMB+1)-Y(:,1)*ones(1,k_LMB);
        [Xr,learned]=rotate(C,Xr,Dz(close,I(2:k_LMB+1)),learning);
        if(~learned)
            nextl=min(l*2,length(index)-1);
            if(nextl>l)
                l=nextl;
                learning=1;
            else
                learning=0;
            end
        else
            l=k_LMB;
            learning=1;
            K=new-ones(1,length(new));
            X(:,J)=Xr(:,K)+X(:,index(front))*ones(1,length(new));
            if mode==1
                Dc(index,J)=dist(X(:,index)',X(:,J));
                Dc(J,index)=dist(X(:,J)',X(:,index));
                if(length(J)~=1)
                    Dc(J,J)=squareform(pdist(X(:,J)'));
                end
            end
            rear=rear+length(new);
            reduced(J)=1;
            front=front+1;
            index=[index,J];
        end
    end
end
end

