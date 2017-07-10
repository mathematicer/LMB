function [ d_hat,local_d ] = local_dimension( D, k, r )
%determine the dimension of the embedded manifold
%k: neighbourhood size
%r: a ratio between 0 and 1
[n,~]=size(D);
local_d=zeros(1,n);
for i=1:n
    [~,I]=sort(D(i,:));
    IDX=I(1:k+1);                                                      
    D_local=D(IDX,IDX);
    W=D_local.*D_local;
    B=eye(k+1)-ones(k+1,k+1)/(k+1);
    X_local_T_X_local=-0.5*B*W*B;
    [~,eigens]=eigs(X_local_T_X_local,k+1);
    eigens=sum(abs(real(eigens)));
    eigensum=sum(eigens);
    threshold=r*eigensum;
    for d=1:k+1
        if(sum(eigens(1:d))>threshold)
            break;
        end
    end    
    local_d(i)=d;
end
d_hat=round(sum(local_d)/n);