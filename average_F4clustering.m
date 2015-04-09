function F_output=average_F4clustering(x,y)
%x, y are the cell arrays
%each cell in x or y represents a class, containing those objects which
%belong to this class

%x is the ground truth
%y is the predicted results

size_x=length(x);
size_y=length(y);

for i=1:size_x
    for j=1:size_y
        
        TP=length(intersect(x{i},y{j}));
        U=union(x{i},y{j});
        
        F_temp(j)=TP/length(U);
        
        
    end

    F(i)=max(F_temp);
end

F_output=mean(F);

