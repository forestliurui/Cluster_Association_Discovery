function F_u=Algo_interior(A,C,F)

[M, k]=size(F);
nIter=100 ;
iteration=1;
while(iteration<nIter)
    B=(C-A.*(F.^2))*ones(k,1)*(ones(k,1)');
    F=((B.^2+4*(A.*C)).^(1/2)-B)./(2*A);
    I=(B(:,1)<0);
    if(sum(I)~=0)
        F(find(I),:)=((diag(F(find(I),:)*ones(k,1)))^(-1))*F(find(I),:);
    end
    iteration=iteration+1;
end
F_u=F;