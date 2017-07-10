function [X,learned] = rotate(C,Xi,D,learning)
%rotate or reflect the local neighbourhood
[~,k]=size(Xi);
[~,l]=size(C);
L=(D.^2-sum(C.^2)'*ones(1,k)-ones(l,1)*sum(Xi.^2))/2;
U=pinv(C')*L*pinv(Xi);
if(abs(abs(det(U))-1)<0.1)
    learned=1;
elseif(l<=4*k && learning==1)
    learned=0;
    X=[];
else
    U=U*(U'*U)^-0.5;
    learned=1;
end
if(learned)
    X=U*Xi;
end
end

