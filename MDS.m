function [X,E] = MDS(D,d)
%Multidimensional scaling
[n,~]=size(D);
W=D.*D;
XTX=-0.5*(eye(n)-ones(n,n)/n)*W*(eye(n)-ones(n,n)/n);
[U,sigma]=eigs(XTX,d);
E=abs(real(sigma));
sigma=sqrt(abs(real(sigma)));
X=sigma*(real(U))';
end