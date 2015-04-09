function MI_normalized=mutual_information_metric(x,y)
%x, y are the cell arrays
%each cell in x or y represents a class, containing those objects which
%belong to this class

%the more similar these two clusterings represented by x and y, the bigger
%MI_normalized

%when MI_normalized = 1, it means that these two clusterings are identical

size_x=length(x);
size_y=length(y);

M_x=0;
for i=1:size_x
    M_x=length(x{i})+M_x;
end
%M_x is the overall number of objects in x


M_y=0;
for i=1:size_y
    M_y=length(y{i})+M_y;
end
%M_y is the overall number of objects in y

if(M_x~=M_y)
    error('mutual_information_metric: the set of objects from two arrays are different');
else
    M=M_x;
end


for i=1:size_x
    for j=1:size_y
        p(i,j)=length(intersect(x{i},y{j}))/M;
    
    end
end

p_x=sum(p,2);
p_y=sum(p,1)';

H_x=0;
for i=1:size_x
    if(p_x(i)~=0)
        H_x=H_x-p_x(i)*log(p_x(i));
    end
end

H_y=0;
for i=1:size_y
    if(p_y(i)~=0)
        H_y=H_y-p_y(i)*log(p_y(i));
    end
end


MI=0;
for i=1:size_x
    for j=1:size_y
        if(p(i,j)==0)
            MI=MI;
        else
            MI=MI+p(i,j)*log(p(i,j)/(p_x(i)*p_y(j)));
        end
    end
end

MI_normalized=MI/max(H_x,H_y);

