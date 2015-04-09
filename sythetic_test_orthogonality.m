

clear;


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
W_i=[1;2;3;4;5;1;2;12;14;13;16];
W_j=[12;13;14;16;11;3;4;7;18;20;19];
W_v=[4;4;4;3;4;1;1;4;4;4;4];

%%synthetic graph 4 (b)
% W_i=[1;2;3;4;5;1;2;12;14;13;16;12;13;14;16];
% W_j=[12;13;14;16;11;3;4;7;18;20;19;1;3;4;2];
% W_v=[4;4;4;3;4;1;1;4;4;4;4;1;1;1;1];

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


intensity_noise=0.3;

W_a_c=sparse(W_i,W_j,W_v,17,20);

% seed = 12345;               % set the seed for random number generator   
% rng(seed);

%add noise to association matrix
W_a_c=W_a_c+max(max(W_a_c))*intensity_noise*rand(size(W_a_c));

% %add noise to association matrix by eliminating nonzero edges
% [r_W,c_W]=find(W_a_c);
% len_W=length(r_W);
% num_rand=randperm(len_W);
% index_remove_noise=num_rand(1:floor(intensity_noise*len_W));
% for i=1:length(index_remove_noise)
%     W_a_c(r_W(index_remove_noise(i)),c_W(index_remove_noise(i)))=0;
% end


k1=4; 
k2=5;  
alpha=1; 
beta=1;
eta1=1;
threshold=0.001;
iteration=1;

% [M0,N0]=size(A0);
% [M1,N1]=size(A1);


%for method_iteration={'NMTFOC','DNMTF','RCC'}
for method_iteration={'NMTFOC','NMTF_BO','DNMTF','RCC'}
method=method_iteration{1,1};
switch method
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

% figure;
% hold on;
% plot(FPR(:,1),TPR(:,1),'--');
% plot(FPR(:,2),TPR(:,2),':');
% plot(FPR(:,3),TPR(:,3),'-');
% legend('NMTFOC','DNMTF','RCC');
% xlabel('FPR');
% ylabel('TPR');
% hold off;


% 
% Theta0=diag(sum(A0,2))-A0;
% Theta1=diag(sum(A1,2))-A1;
% 
% H0=10*rand(M0, k0);
% H1=10*rand(M1, k1);
% S=10*rand(k0,k1);
% % S(1,1)=30;
% % S(2,2)=30;
% W_p=(abs(W)+W)./2;
% W_n=(abs(W)-W)./2;
% Theta0_n=(abs(Theta0)-Theta0)./2;
% Theta0_p=(abs(Theta0)+Theta0)./2;
% Theta1_n=(abs(Theta1)-Theta1)./2;
% Theta1_p=(abs(Theta1)+Theta1)./2;
% 
% S_u=S.*((H0'*W_p*H1)./((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )).^(1/2);
% H0_u=H0.*((W_p*H1*S'+alpha*Theta0_n*H0)./((H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0)))).^(1/2);
% H1_u=H1.*((W_p'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1) )).^(1/2);
% 
% % S_n= (((H0'*W*H1)-0.5*eta1)./(H0'*H0*S*(H1')*H1)).*S;
% % H0_n= (W*H1*(S'))./(alpha*(Theta0*H0)+(H0*S*(H1')*H1*S'  ) ).*H0;
% % H1_n= ((W')*H0*S)./(beta*(Theta1*H1)+(H1*(S')*(H0')*H0*S  ) ).*H1;
% 
% S_res=norm(S_u-S);
% H0_res=norm(H0_u-H0);
% H1_res=norm(H1_u-H1);
% 
% while((S_res>threshold)||(H0_res>threshold)||(H1_res>threshold))
%     iteration=iteration+1;
%     S=S_u;
%     H0=H0_u;
%     H1=H1_u;
% %     S_n= (((H0'*W*H1)-0.5*eta1)./((H0')*H0*S*(H1')*H1)).*S;
% %     H0_n= ((W*H1*(S'))./(alpha*(Theta0*H0)+(H0*S*(H1')*H1*(S')  ) )).*H0;
% %     
% %     H0_n= H0_n./(sum(H0_n,2)*ones(1,k0));
% %     
% %     H1_n= (((W')*H0*S)./(beta*(Theta1*H1)+(H1*(S')*(H0')*H0*S  ) )).*H1;
% %     
% %     H1_n= H1_n./(sum(H1_n,2)*ones(1,k1));
% 
%     S_u=S.*((H0'*W_p*H1)./((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )).^(1/2);
%     
% %    S_u=S_u/sum(sum(S_u));
%      
%    H0_u=H0.*((W_p*H1*S'+alpha*Theta0_n*H0)./((H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0)))).^(1/2);
%     
%     H0_u= H0_u./(sum(H0_u,2)*ones(1,k0));
%     
%     H1_u=H1.*((W_p'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1) )).^(1/2);
%     
%     H1_u= H1_u./(sum(H1_u,2)*ones(1,k1));
%     
%     S_res=norm(S_u-S);
%     H0_res=norm(H0_u-H0);
%     
%     H1_res=norm(H1_u-H1);
% end
% 
% [val_H0, ind_H0]=max(H0,[],2);
% [val_H1, ind_H1]=max(H1,[],2);
% 
% f_res=norm(W-(H0)*S*(H1'))+alpha*trace((H0')*Theta0*H0)+beta*trace((H1')*Theta1*H1)+eta1*norm(S(:),1)
