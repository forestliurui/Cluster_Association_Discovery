function output=voting_strategy(input)
%find the number appearing most frequently in the vector of input

N=length(input);

for i=1:max(input)
    count_v(i)=0;
end

for i=1:N
    count_v(input(i))=count_v(input(i))+1;
end
[num output]=max(count_v);