%test on mutual information
clear;
G1_ground_truth={[1,2,3,4],[12,13,14,16],[9,11],[5,6,8,7,10],[15,17]};
G1_class_predicted_NMTFOC={[12,13,14,15,16,17],[1,2,3,4],[5,6,7,8,10],[9,11]};
G2_class_predicted_NMTFOC={[18,19,20],[12,13,14,15,16,17],[9,11],[]};

G1_class_predicted_NMTF_Chris={[2,4,13],[5,6,7,8,9,10,11,14,15,17],[1,3,12],[16]};
G1_class_predicted_DNMTF={[1:17]};
G1_class_predicted_RCC={[1,2,6,9,10,16],[3,5,7,8],[4,15],[11,12,13,14,17]};



MI=mutual_information_metric(G1_class_predicted_NMTFOC,G1_ground_truth);
MI_chris=mutual_information_metric(G1_class_predicted_NMTF_Chris,G1_ground_truth);
MI_DNMTF=mutual_information_metric(G1_class_predicted_DNMTF,G1_ground_truth);
MI_RCC=mutual_information_metric(G1_class_predicted_RCC,G1_ground_truth);


%%

a={[1]}



%%
clear;
p=[4,0,2;
    0,5,2;
    2,2,4];
p=p/sum(sum(p));
[size_x, size_y]=size(p);

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