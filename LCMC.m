function [ ave ] = LCMC( k,D,X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,npoints]=size(X);
INTESECT=zeros(1,npoints);
for i=1:npoints
    [~,IDX]=sort(D(i,:));
    I1=IDX(2:k+1);
    IDX=knnsearch(X',X(:,i)','k',k+1);
    I2=IDX(2:k+1);
    INTESECT(i)=length(intersect(I1,I2));
end
ave=mean(INTESECT)/k;
end