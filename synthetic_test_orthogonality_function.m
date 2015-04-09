function result_4method=synthetic_test_orthogonality_function(intensity_noise)




%% Beginning-- the input for graph without feature 
A1_i=[1;1;1;2;3;4;5;6;5;7;4;7;9;5;9;11;12;12;13;14;13;14;13];
A1_j=[2;3;4;4;4;5;6;7;8;8;9;10;11;12;12;13;13;14;14;15;16;16;17];
A1_v=[4;4;3;3;4;1;1;1;1;1;1;1;4;1.5;1;1;4;4;3;1;3.5;2;2];

A1_i_d=(1:17)';
A1_j_d=(1:17)';
A1_v_d=3.5*ones(17,1);




A2_i=[1;1;1;2;3;4;5;6;5;7;4;7;9;5;9;11;12;12;13;14;13;14;13;7;7;18;18;19];
A2_j=[2;3;4;4;4;5;6;7;8;8;9;10;11;12;12;13;13;14;14;15;16;16;17;18;20;20;19;20];
A2_v=[4;4;3;3;4;1;1;1;1;1;1;1;4;1.5;1;1;4;4;3;1;3.5;2;2;4;4;4;3;3.5];

A2_i_d=(1:20)';
A2_j_d=(1:20)';
A2_v_d=3.5*ones(20,1);

% A_i=[A_i;A_i_d];
% A_j=[A_j;A_j_d];
% A_v=[A_v;A_v_d];
% A1_i=[A1_i;A1_i_d];
% A1_j=[A1_j;A1_j_d];
% A1_v=[A1_v;A1_v_d];

A1=sparse(A1_i,A1_j,A1_v,17,17);
A2=sparse(A1_i,A1_j,A1_v,20,20);
A1=A1+A1';
A2=A2+A2';

%%1
% W_i=[1;2;3;4;5;1;2];
% W_j=[12;13;14;16;11;3;4];
% W_v=[4;4;4;3;4;1;1];

%%2
% W_i=[1;2;3;4;5;1;2;12;14;13;16];
% W_j=[12;13;14;16;11;3;4;7;18;20;19];
% W_v=[4;4;4;3;4;1;1;1;1;1;1];

%%3 (a)
% W_i=[1;2;3;4;5;1;2;12;14;13;16];
% W_j=[12;13;14;16;11;3;4;7;18;20;19];
% W_v=[4;4;4;3;4;1;1;4;4;4;4];

%%synthetic graph 4 (b)
W_i=[1;2;3;4;5;1;2;12;14;13;16;12;13;14;16];
W_j=[12;13;14;16;11;3;4;7;18;20;19;1;3;4;2];
W_v=[4;4;4;3;4;1;1;4;4;4;4;1;1;1;1];

%%synthetic graph 4 (c)
% W_i=[1;2;3;4;5;1;2;12;14;13;16;12;13;14;16;  10;6;7;7];
% W_j=[12;13;14;16;11;3;4;7;18;20;19;1;3;4;2; 10;9;1;17];
% W_v=[4;4;4;3;4;1;1;4;4;4;4;1;1;1;1;   1;1;1;1];

% C_association_true=zeros(17,20);
% C_association_predicted=zeros(17,20);
% 
% C_association_true([1,2,3,4],[12,13,14,16])=1;
% C_association_true([12,13,14,16],[7,18,19,20])=1;
% C_association_true([1,2,3,4],[1,2,3,4])=1;
% C_association_true([5],[11])=1;

domain1_truth{1}=[1;2;3;4];
domain2_truth{1}=[12;13;14;16];
domain1_truth{2}=[12;13;14;16];
domain2_truth{2}=[7;18;19;20];




% W_i=[1;2;3;4;5;1;2;12;14;13;16;12;13;14;16];
% W_j=[12;13;14;16;11;3;4;7;18;20;19;1;3;4;2];
% W_v=[4;4;4;3;4;1;1;4;4;4;4;4;4;4;4];


% intensity_noise=0.3;

W_a_c=sparse(W_i,W_j,W_v,17,20);

% seed = 12345;               % set the seed for random number generator   
% rng(seed);

%add noise to association matrix
W_a_c=W_a_c+max(max(W_a_c))*intensity_noise*rand(size(W_a_c));
%a different way to add noise to association matrix--flip the elements
% [M_Wac,N_Wac]=size(W_a_c);
% rand_MNac=randperm(M_Wac*N_Wac);
% 
%  W_a_c(rand_MNac(1:floor(M_Wac*N_Wac*intensity_noise)))=4*(~W_a_c(rand_MNac(1:floor(M_Wac*N_Wac*intensity_noise))));%flip intensity_noise percentage of elements in W


% %add noise to association matrix by eliminating nonzero edges
% [r_W,c_W]=find(W_a_c);
% len_W=length(r_W);
% num_rand=randperm(len_W);
% index_remove_noise=num_rand(1:floor(intensity_noise*len_W));
% for i=1:length(index_remove_noise)
%     W_a_c(r_W(index_remove_noise(i)),c_W(index_remove_noise(i)))=0;
% end

%% END-- the input for graph without feature 

%% Beginning--the input for graph with feature
F1=csvread('sythetic_feature\\synthetic_A1.csv');


W=5*csvread('sythetic_feature\\synthetic_W.csv');
[M_W,N_W]=size(W);
rand_MN=randperm(M_W*N_W);

 W(rand_MN(1:floor(M_W*N_W*intensity_noise)))=5*(~W(rand_MN(1:floor(M_W*N_W*intensity_noise))));%flip intensity_noise percentage of elements in W
%     k=intensity_noise;
%    W=W+intensity_noise*rand(size(W))*max(max(W));
F2=F1*W;
%    F2=F1*W+max(max(F1*W))*intensity_noise*(rand(size(F1*W))-0);
% SNR(k)=max(max(F1*W))/intensity_noise;

F1=F1;
for i=1:size(F1,2)
    A1_ft(:,i)=sum((F1(:,i)*ones(1, size(F1,2))).*F1,1);
end

for i=1:size(F2,2)
    A2_ft(:,i)=sum((F2(:,i)*ones(1, size(F2,2))).*F2,1);
end
% 
% domain1_truth{1}=[1;2;3;4];
% domain2_truth{1}=[1;2;3;4];
% domain1_truth{2}=[5;6;7;8];
% domain2_truth{2}=[5;6;7;8];
% domain1_truth{3}=[9;10;11;12];
% domain2_truth{3}=[9;10;11;12];
%% END-- the input for graph with feature

k1=3; 
k2=3;  
alpha=1; 
beta=1;
eta1=1;
threshold=0.001;
iteration=1;

% [M0,N0]=size(A0);
% [M1,N1]=size(A1);


for method_iteration={'NMTFOC','NMTF_BO','DNMTF','RCC'}
% for method_iteration={'NMTFOC_ft'}
method=method_iteration{1,1};
switch method
    case 'NMTFOC_ft'
        [ind_H1, ind_H2, S, H1, H2, W]=NMTFOC_ft(A1_ft,A2_ft,F1,F2,k1,k2);
        
    case 'DRCC',
        [ind_H1, ind_H2, S, H1, H2]=DRCC(W_a_c, A1, A2, k1, k2);
    case 'DNMTF'
        [ind_H1, ind_H2, S, H1, H2]=DNMTF(W_a_c, A1, A2, k1,k2);
    case 'NMTFOC' 
        [ind_H1, ind_H2, S, H1, H2]=NMTFOC(W_a_c, A1, A2, k1,k2);
        
        case 'NMTF_BO'
        [ind_H1, ind_H2, S, H1, H2]=NMTF_BO(W_a_c, A1, A2, k1,k2);
        
    case 'RCC'
        [H1,H2,S]=RCC(W_a_c, A1, A2, k1, k2);
        [val_H1, ind_H1]=max(H1,[],2);
        [val_H2, ind_H2]=max(H2,[],2);   
    case 'NMTFMG'
%         F_cell{1,1}=F1;
%         F_cell{2,2}=F2;
%         A_cell{1,1}=A0/mean(mean(A0));
%         A_cell{2,2}=A1/mean(mean(A1));
%         c=[10,10];
        display('enter NMTFMG\n');
 %       [ind_H,  S_rec, H_rec, W_rec]=NMTFMG( A_cell, F_cell,c, W_cell);
         [ind_H,  S_rec, H_rec, W_rec]=NMTFMG_noregression( A_cell, F_cell,c, W_cell);
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
if(strcmp(method,'NMTFMG'))
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

timestamp=[int2str(year(now)),'_',int2str(month(now)),'_',int2str(day(now)),'_',int2str(hour(now)),'_',int2str(minute(now))];
% save(['results_DBLP_labeled_',method,'_',timestamp,'.mat'],'ind_H','ind_H1','ind_H2','S_rec','S','H_rec','H1','H2','W_rec','W');

% if(strcmp(method,'NMTFMG'))
%     save(['results_syn_',method,'_',timestamp,'.mat'],'ind_H','ind_H1','ind_H2','S_rec','S','H_rec','H1','H2','W_rec','W');
% else
%     save(['results_syn_',method,'_',timestamp,'.mat'],'ind_H1','ind_H2','S','H1','H2');
% end

% filename=['results_syn_',method,'_',timestamp,'.txt'];
% save(['results_yeast_',method,'_',timestamp,'.mat'],'ind_H','ind_H1','ind_H2','S_rec','S','H_rec','H1','H2','W_rec','W');
% 
% filename=['results_yeast_',method,'_',timestamp,'.txt'];

% ffile=fopen(filename,'wt');

C_association_predicted=zeros(17,20);
S_temp=S;
for i=1:2
[rows, cols]=find(max(max(S_temp))==S_temp);
rows=rows(1);
cols=cols(1);

domain1_predicted{i}=find(ind_H1==rows);
domain2_predicted{i}=find(ind_H2==cols);

display(['The ', num2str(i),'-th largest pair is Row: ', num2str(rows),' Col: ',num2str(cols)]);

% C_association_predicted(find(ind_H1==rows),find(ind_H2==cols))=1;
% fprintf(ffile,'%s\n',['right cluster ',num2str(i)]);
% fprintf(ffile,'%d\n',find(ind_H2==cols));
% fprintf(ffile,'\n\n\n');
% fprintf(ffile,'%s\n',['left cluster ',num2str(i)]);
% fprintf(ffile,'%d\n',find(ind_H1==rows)   );
% fprintf(ffile,'\n');
A1_average=tri_dense_eval(full(A1),full(A2),W,find(ind_H1==rows),[]);
A2_average=tri_dense_eval(full(A1),full(A2),W,[],find(ind_H2==cols));
W_average=sum(sum(W(find(ind_H1==rows),find(ind_H2==cols))))/(sum(ind_H1==rows)+sum(ind_H2==cols));
% 
% fprintf(ffile,'the majority of authors from author cluster %d are from %d\n',i,voting_strategy(authorlist_label(find(ind_H1==rows))+1)-1);
% 
% fprintf(ffile,'the percent of persons from database: %f; Data Mining: %f; AI: %f; IR: %f\n', sum(authorlist_label(find(ind_H1==rows))==0)/length(authorlist_label(find(ind_H1==rows))),sum(authorlist_label(find(ind_H1==rows))==1)/length(authorlist_label(find(ind_H1==rows))), sum(authorlist_label(find(ind_H1==rows))==2)/length(authorlist_label(find(ind_H1==rows))), sum(authorlist_label(find(ind_H1==rows))==3)/length(authorlist_label(find(ind_H1==rows))) );
% display(['the percent of persons from database: ',num2str(sum(authorlist_label(find(ind_H1==rows))==0)/length(authorlist_label(find(ind_H1==rows)))),'; Data Mining: ',num2str(sum(authorlist_label(find(ind_H1==rows))==1)/length(authorlist_label(find(ind_H1==rows)))),'; AI: ',num2str(sum(authorlist_label(find(ind_H1==rows))==2)/length(authorlist_label(find(ind_H1==rows)))),'; IR: ',num2str(sum(authorlist_label(find(ind_H1==rows))==3)/length(authorlist_label(find(ind_H1==rows)))),'\n'] );
% 
% fprintf(ffile,'the percent(multiple label) of persons from database: %f; Data Mining: %f; AI: %f; IR: %f\n', sum(multi_label(find(ind_H1==rows),1)==1)/length(authorlist_label(find(ind_H1==rows))),sum(multi_label(find(ind_H1==rows),2)==1)/length(authorlist_label(find(ind_H1==rows))),sum(multi_label(find(ind_H1==rows),3)==1)/length(authorlist_label(find(ind_H1==rows))), sum(multi_label(find(ind_H1==rows),4)==1)/length(authorlist_label(find(ind_H1==rows))) );
% display(['the percent(multiple label) of persons from database: ',num2str(sum(multi_label(find(ind_H1==rows),1)==1)/length(authorlist_label(find(ind_H1==rows)))),'; Data Mining: ',num2str(sum(multi_label(find(ind_H1==rows),2)==1)/length(authorlist_label(find(ind_H1==rows)))),'; AI: ',num2str(sum(multi_label(find(ind_H1==rows),3)==1)/length(authorlist_label(find(ind_H1==rows)))),'; IR: ',num2str(sum(multi_label(find(ind_H1==rows),4)==1)/length(authorlist_label(find(ind_H1==rows)))),'\n'] );


% fprintf(ffile, 'the tri-dense score for A1 is %f, for A2 is %f, for W is %f, the overall average is %f. The previous way to calculation get %f\n\n\n',A1_average,A2_average,W_average,tri_dense_eval(full(A1),full(A2),W,find(ind_H1==rows),find(ind_H2==cols)));
% fprintf(ffile,'\n\n');
S_temp(rows,cols)=0;
end
% 
% asso_measure=matrix_consistency4association(C_association_true,C_association_predicted);
% fprintf(ffile, 'measure of association: %f\n', asso_measure);
% display(['measure of association: ',num2str(asso_measure)]);


% fclose(ffile);



F_score=clustering_association_measure(domain1_predicted, domain2_predicted, domain1_truth,domain2_truth);

% fprintf(ffile, 'F_score: %f\n', F_score);
fprintf(1, 'F_score: %f\n', F_score);

% fclose(ffile);

result_4method(iteration)=F_score;

% 
% [TPR(:,iteration), FPR(:,iteration), AUC(iteration)]=roc_curve(full(W),full(W_a_c)>0);
% if(strcmp(method,'NMTFOC'))
%     AUC_NMTFOC(k)=AUC;
% elseif(strcmp(method,'RCC'))
%     AUC_RCC(k)=AUC;
% else
%     AUC_DNMTF(k)=AUC;
% end

iteration=iteration+1;
end

