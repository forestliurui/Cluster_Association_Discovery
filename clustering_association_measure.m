function result_measure=clustering_association_measure(C1_predicted, C2_predicted, C1_true, C2_true)

if(length(C1_predicted)~=length(C2_predicted))
    error('clustering_association_measure: predicted data does not form pairs');
end

if(length(C1_true)~=length(C2_true))
    error('clustering_association_measure: groundtruth data does not form pairs');
end

max_C1_true=-inf;

for i=1:length(C1_true)
    max_C1_true=max(max_C1_true, max(C1_true{i}));

end


for i=1:length(C2_predicted)
    C2_predicted{i}=C2_predicted{i}+max_C1_true+10;
end
for i=1:length(C2_true)
    C2_true{i}=C2_true{i}+max_C1_true+10;
end

for i=1:length(C1_predicted)
    C_predicted{i}=union(C1_predicted{i}, C2_predicted{i});
    if(size(C_predicted{i},2)~=1)
        C_predicted{i}=C_predicted{i}';
    end
end

for i=1:length(C1_true)
    C_true{i}=union(C1_true{i}, C2_true{i});
end

all_element_C=[];
for i=1:length(C_predicted)
    all_element_C=[all_element_C;C_predicted{i}];
end

all_element_C=reshape(all_element_C, size(all_element_C,1)*size(all_element_C,2),1);

i=1;
for j=1:length(C_true)
    C_true4MI_temp=intersect(C_true{j},all_element_C);
    if(~isempty(C_true4MI_temp))
        C_true4MI{i}=C_true4MI_temp;
        i=i+1;
    end
        
end

% MI_normalized=mutual_information_metric(C_predicted,C_true4MI); %mutual
% infomration
F_output=average_F4clustering(C_true,C_predicted);  %CAM

result_measure=F_output;