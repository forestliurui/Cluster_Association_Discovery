function MI=mutual_information_matrix(W,H1,H2)

%normalize W
W=W/sum(sum(W));

num_clus1=max(H1);
num_clus2=max(H2);

for i=1:num_clus1
    for j=1:num_clus2
        p(i,j)= sum(sum(W(find(H1==i),find(H2==j))));
    end
end

p_x=sum(p,2);
p_y=sum(p,1)';


MI=0;
for i=1:num_clus1
    for j=1:num_clus2
        if(p(i,j)==0)
            MI=MI;
        else
            MI=MI+p(i,j)*log(p(i,j)/(p_x(i)*p_y(j)));
        end
    end
end

MI_normalized=MI;

