function C=gdivide_revised(A,B)

if((size(A,1)~=size(B,1))||(size(A,2)~=size(B,2)))
    error('gdivide_revised: the dimensions of two input matrix don\''t match');
end

[M,N]=size(B);

for i=1:M
    for j=1:N
        if(B(i,j)==0)
            C(i,j)=0;
        else
            C(i,j)=A(i,j)/B(i,j);
        end
    end
end

    