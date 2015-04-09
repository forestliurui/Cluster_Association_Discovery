clear;
% load('DBLP_labeled_feature.mat');
% load('DBLP_labeled_feature_removeDMandSTOP.mat');
load('DBLP_labeled_feature_removeandSTOP.mat');
% load('DBLP_labeled_feature_removeIRandSTOP.mat');
% load('DBLP_labeled_feature_removeDMIRandSTOP.mat');
% load('DBLP_labeled_feature_stopwords_removed.mat');

for i=1:size(A0_a_a,1)
    A0_a_a(i,i)=0;
end

author_remove_caused_by_label=[];
author_remain_caused_by_label=1:size(A0_a_a,1);

%get the index of authos to be removed with label DM and IR 
% author_remove_caused_by_label=find((authorlist_label==1)|(authorlist_label==3));
% author_remain_caused_by_label=find(~((authorlist_label==1)|(authorlist_label==3)));

% author_remove_caused_by_label=find((authorlist_label==1)); %remove DM
% author_remain_caused_by_label=find(~((authorlist_label==1)));

%get the index of authors who has no co-authorship with other if 3 papers are used
%as the minimal threshold for co-authorship
author_remain_ind=intersect(author_remain_caused_by_label,find(sum(A0_a_a>3,2)));
author_remove_ind=union(author_remove_caused_by_label,find(~sum(A0_a_a>3,2)));

%remove those authoers from our input dataset
authorlist_name(author_remove_ind)=[];
authorlist_label=authorlist_label(author_remain_ind);
authorlist_ID=authorlist_ID(author_remain_ind);
A0_a_a=A0_a_a(author_remain_ind,author_remain_ind);
F_t_a=F_t_a(:,author_remain_ind);
W_a_c=W_a_c(author_remain_ind,:); 

intensity_noise=0.1;
% intensity_noise=0.15;
%add noise to association matrix
% W_a_c=W_a_c+max(max(W_a_c))*0.2*rand(size(W_a_c));

% seed = 12345;               % set the seed for random number generator   
% rng(seed);
%add noise to association matrix by eliminating nonzero edges
[r_W,c_W]=find(W_a_c);
len_W=length(r_W);
num_rand=randperm(len_W);
index_remove_noise=num_rand(1:floor(intensity_noise*len_W));
for i=1:length(index_remove_noise)
    W_a_c(r_W(index_remove_noise(i)),c_W(index_remove_noise(i)))=0;
end

W_index_remove_author=find(~sum(W_a_c,2)); %sum along all conferences
W_index_remove_conf=find(~sum(W_a_c,1)); %sum along all authors

W_a_c(W_index_remove_author,:)=[];
W_a_c(:,W_index_remove_conf)=[];

authorlist_name(W_index_remove_author)=[];
authorlist_label(W_index_remove_author)=[];
authorlist_ID(W_index_remove_author)=[];
A0_a_a(W_index_remove_author,:)=[];
A0_a_a(:,W_index_remove_author)=[];
F_t_a(:,W_index_remove_author)=[];



multi_label_threshold=3;
%multi_label(:,1) denotes the label for data base,  multi_label(:,2) denotes the label for data mining
for i=1:4
    multi_label_temp=sum(W_a_c(:,find(conf_label==i-1)),2);
    multi_label(:,i)=(full(multi_label_temp)>multi_label_threshold);
end

for i=1:size(multi_label,1)
    multi_label_test(i)=multi_label(i,authorlist_label(i)+1);
end


W_cell{1,2}=W_a_c/mean(mean(W_a_c));
W_cell{2,1}=W_a_c'/mean(mean(W_a_c));
F_cell{1,1}=F_t_a/mean(mean(F_t_a));
F_cell{2,2}=F_t_c/mean(mean(F_t_c));


% for i=1:size(F_ta,2)
%     A1(:,i)=sum((F1(:,i)*ones(1, size(F1,2))).*F1,1);
% end

A1=A0_a_a;

[conf_l_ID conf_l_label conf_l_list]=textread('DBLP_four_area_labeled\\conf_label.txt','%d%d%s','delimiter','\t');

% A2_label=zeros(length(conf_label));
% for i=0:length(conf_l_label)
%     A2_label(conf_label==i,conf_label==i)=5;
% end
% balance_label=10;

for i=1:size(F_t_c,2)
    A2(:,i)=sum((F_t_c(:,i)*ones(1, size(F_t_c,2))).*F_t_c,1);
end



A_cell{1,1}=A1/mean(mean(A1));
% A_cell{2,2}=A2/mean(mean(A2))+balance_label*A2_label;
A_cell{2,2}=A2/mean(mean(A2));
% A_cell{2,2}=A2_label;
c=[10,2];
k2=4;
k1=4;

for method_iteration={'NMTFOC','NMTF_BO','DNMTF','RCC'}
method=method_iteration{1,1};
switch method
    case 'DRCC',
        [ind_H1, ind_H2, S, H1, H2]=DRCC(W_a_c, A1, A2, k1, k2);
    case 'DNMTF'
        [ind_H1, ind_H2, S, H1, H2]=DNMTF(W_a_c, A1, A2, k1,k2);
    case 'NMTF_BO'
        [ind_H1, ind_H2, S, H1, H2]=NMTF_BO(W_a_c, A1, A2, k1,k2);
        
    case 'NMTFOC' 
        [ind_H1, ind_H2, S, H1, H2]=NMTFOC(W_a_c, A1, A2, k1,k2);
    case 'NMTFOC_ft' 
        [ind_H1, ind_H2, S, H1, H2,W]=NMTFOC_ft(A1, A2,F_cell{1,1},F_cell{2,2},k1,k2,W_a_c);
    case 'RCC'
        [H1,H2,S]=RCC(W_a_c, A1, A2, k1, k2);
        [val_H1, ind_H1]=max(H1,[],2);
        [val_H2, ind_H2]=max(H2,[],2);   
    case 'NMTFMG_p'
%         F_cell{1,1}=F1;
%         F_cell{2,2}=F2;
%         A_cell{1,1}=A0/mean(mean(A0));
%         A_cell{2,2}=A1/mean(mean(A1));
%         c=[10,10];
        display('enter NMTFMG\n');
        [ind_H,  S_rec, H_rec, W_rec]=NMTFMG_p( A_cell, F_cell,c, W_cell);
%         [ind_H,  S_rec, H_rec, W_rec]=NMTFMG( A_cell, F_cell,c);
        display('get out of NMTFMG');
    otherwise,
        ;
end


% method='NMTFMG';
% display('enter NMTFMG\n');
% % [ind_H,  S_rec, H_rec, W_rec]=NMTFMG( A_cell, F_cell,c,W_cell);
% [ind_H,  S_rec, H_rec, W_rec]=NMTFMG( A_cell, F_cell,c);
% display('get out of NMTFMG');
if(strcmp(method,'NMTFMG_p'))
S=S_rec(1:c(1),c(1)+1:c(1)+c(2));
H1=H_rec(1:size(W_cell{1,2},1),1:c(1));
H2=H_rec(size(W_cell{1,2},1)+1:size(W_cell{1,2},1)+size(W_cell{1,2},2),c(1)+1:c(1)+c(2));
[val_H1, ind_H1]=max(H1,[],2);
[val_H2, ind_H2]=max(H2,[],2);
W=full(W_rec(1:size(W_cell{1,2},1),size(W_cell{1,2},1)+1:size(W_cell{1,2},1)+size(W_cell{1,2},2)));
else
    W=H1*S*(H2');
end
%%

for i=0:3
    conf_class_truth4MI{i+1}=find(conf_label==i);
    conf_class_predicted4MI{i+1}=find(ind_H2==(i+1));
    
    author_class_truth4MI{i+1}=find(authorlist_label==i);
    author_class_predicted4MI{i+1}=find(ind_H1==(i+1));     
end

MI_conf=mutual_information_metric(conf_class_predicted4MI,conf_class_truth4MI);
MI_author=mutual_information_metric(author_class_predicted4MI,author_class_truth4MI);

fprintf(1,'MI of W: %f\n',mutual_information_matrix(W_a_c,ind_H1,ind_H2));

timestamp=[int2str(year(now)),'_',int2str(month(now)),'_',int2str(day(now)),'_',int2str(hour(now)),'_',int2str(minute(now))];
% save(['results_DBLP_labeled_',method,'_',timestamp,'.mat'],'ind_H','ind_H1','ind_H2','S_rec','S','H_rec','H1','H2','W_rec','W');

if(strcmp(method,'NMTFMG_p'))
    save(['results_DBLP_labeled_',method,'_',timestamp,'.mat'],'ind_H','ind_H1','ind_H2','S_rec','S','H_rec','H1','H2','W_rec','W');
elseif(strcmp(method,'NMTFOC_ft'))
    save(['results_DBLP_labeled_',method,'_',timestamp,'.mat'],'ind_H1','ind_H2','S','H1','H2','W');
else
    save(['results_DBLP_labeled_',method,'_',timestamp,'.mat'],'ind_H1','ind_H2','S','H1','H2');
end


filename=['output_DBLP_four_area\DBLP_noise_MI_value\results_DBLP_labeled_',method,'_',timestamp,'.txt'];
% save(['results_yeast_',method,'_',timestamp,'.mat'],'ind_H','ind_H1','ind_H2','S_rec','S','H_rec','H1','H2','W_rec','W');
% 
% filename=['results_yeast_',method,'_',timestamp,'.txt'];

ffile=fopen(filename,'wt');


S_temp=S;
for i=1:4
[rows, cols]=find(max(max(S_temp))==S_temp);
display(['The ', num2str(i),'-th largest pair is Row: ', num2str(rows),' Col: ',num2str(cols)]);

clustering_association_predicted_author{i}= find(ind_H1==rows) ;

clustering_association_predicted_conf{i}= find(ind_H2==cols) ;

% clustering_association_true_author{i}= find(ind_H1==rows) ;
% 
% clustering_association_true_conf{i}= find(ind_H2==cols) ;


fprintf(ffile, 'MI_conf: %f; MI_author: %f\n', MI_conf,MI_author);
fprintf(ffile,'%s\n',['conference cluster ',num2str(i)]);
fprintf(ffile,'%s\n',conflist{find(ind_H2==cols)});
fprintf(ffile,'\n\n\n');
fprintf(ffile,'%s\n',['author cluster ',num2str(i)]);
fprintf(ffile,'%s\n',authorlist_name{  find(ind_H1==rows)  } );
fprintf(ffile,'\n');
A1_average=tri_dense_eval(A_cell{1,1},A_cell{2,2},W,find(ind_H1==rows),[]);
A2_average=tri_dense_eval(A_cell{1,1},A_cell{2,2},W,[],find(ind_H2==cols));
W_average=sum(sum(W(find(ind_H1==rows),find(ind_H2==cols))))/(sum(ind_H1==rows)+sum(ind_H2==cols));

% fprintf(ffile,'the majority of authors from author cluster %d are from %d\n',i,voting_strategy(authorlist_label(find(ind_H1==rows))+1)-1);

fprintf(ffile,'the percent of persons from database: %f; Data Mining: %f; AI: %f; IR: %f\n', sum(authorlist_label(find(ind_H1==rows))==0)/length(authorlist_label(find(ind_H1==rows))),sum(authorlist_label(find(ind_H1==rows))==1)/length(authorlist_label(find(ind_H1==rows))), sum(authorlist_label(find(ind_H1==rows))==2)/length(authorlist_label(find(ind_H1==rows))), sum(authorlist_label(find(ind_H1==rows))==3)/length(authorlist_label(find(ind_H1==rows))) );
display(['the percent of persons from database: ',num2str(sum(authorlist_label(find(ind_H1==rows))==0)/length(authorlist_label(find(ind_H1==rows)))),'; Data Mining: ',num2str(sum(authorlist_label(find(ind_H1==rows))==1)/length(authorlist_label(find(ind_H1==rows)))),'; AI: ',num2str(sum(authorlist_label(find(ind_H1==rows))==2)/length(authorlist_label(find(ind_H1==rows)))),'; IR: ',num2str(sum(authorlist_label(find(ind_H1==rows))==3)/length(authorlist_label(find(ind_H1==rows)))),'\n'] );



% fprintf(ffile,'the percent(multiple label) of persons from database: %f; Data Mining: %f; AI: %f; IR: %f\n', sum(multi_label(find(ind_H1==rows),1)==1)/length(authorlist_label(find(ind_H1==rows))),sum(multi_label(find(ind_H1==rows),2)==1)/length(authorlist_label(find(ind_H1==rows))),sum(multi_label(find(ind_H1==rows),3)==1)/length(authorlist_label(find(ind_H1==rows))), sum(multi_label(find(ind_H1==rows),4)==1)/length(authorlist_label(find(ind_H1==rows))) );
% display(['the percent(multiple label) of persons from database: ',num2str(sum(multi_label(find(ind_H1==rows),1)==1)/length(authorlist_label(find(ind_H1==rows)))),'; Data Mining: ',num2str(sum(multi_label(find(ind_H1==rows),2)==1)/length(authorlist_label(find(ind_H1==rows)))),'; AI: ',num2str(sum(multi_label(find(ind_H1==rows),3)==1)/length(authorlist_label(find(ind_H1==rows)))),'; IR: ',num2str(sum(multi_label(find(ind_H1==rows),4)==1)/length(authorlist_label(find(ind_H1==rows)))),'\n'] );


% fprintf(ffile, 'the tri-dense score for A1 is %f, for A2 is %f, for W is %f, the overall average is %f. The previous way to calculation get %f\n\n\n',A1_average,A2_average,W_average,tri_dense_eval(full(A_cell{1,1}),full(A_cell{2,2}),W,find(ind_H1==rows),find(ind_H2==cols)));
% fprintf(ffile,'\n\n');
S_temp(rows,cols)=0;
end


% MI_association_normalized=clustering_association_measure(clustering_association_predicted_author, clustering_association_predicted_conf, author_class_truth4MI,conf_class_truth4MI);
F_score=clustering_association_measure(clustering_association_predicted_author, clustering_association_predicted_conf, author_class_truth4MI,conf_class_truth4MI);

fprintf(ffile, 'F_score: %f\n', F_score);
fprintf(1, 'F_score: %f\n', F_score);

fclose(ffile);
end