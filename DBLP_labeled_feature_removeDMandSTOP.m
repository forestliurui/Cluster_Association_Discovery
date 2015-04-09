%%
clear;


[authorlist_ID_input authorlist_label_input authorlist_name_input]=textread('DBLP_four_area_labeled\\author_label.txt','%d%d%s','delimiter','\t');
[p_a1_input p_a2_input]=textread('DBLP_four_area_labeled\\paper_author.txt','%d%d','delimiter','\t');
[p_t1_input p_t2_input]=textread('DBLP_four_area_labeled\\paper_term.txt','%d%d','delimiter','\t');
[termlist_ID_input termlist_input]=textread('DBLP_four_area_labeled\\term.txt','%d%s','delimiter','\t');
[conf_ID_input, conf_label_input, conflist_input]=textread('DBLP_four_area_labeled\\conf_label.txt','%d%d%s','delimiter','\t');
[p_c1_input p_c2_input]=textread('DBLP_four_area_labeled\\paper_conf.txt','%d%d','delimiter','\t');

conf_ID_DB=conf_ID_input(conf_label_input==0);
conf_label_DB=0*ones(length(conf_ID_DB),1);
conflist_DB=conflist_input(conf_label_input==0);
p_c2_DB=p_c2_input(ismember(p_c2_input,conf_ID_DB));
p_c1_DB=p_c1_input(ismember(p_c2_input,conf_ID_DB));
p_t1_DB=p_t1_input(ismember(p_t1_input,p_c1_DB));
p_t2_DB=p_t2_input(ismember(p_t1_input,p_c1_DB));
p_a1_DB=p_a1_input(ismember(p_a1_input,p_c1_DB));
p_a2_DB=p_a2_input(ismember(p_a1_input,p_c1_DB));

conf_ID_DM=conf_ID_input(conf_label_input==1);
conf_label_DM=1*ones(length(conf_ID_DM),1);
conflist_DM=conflist_input(conf_label_input==1);
p_c2_DM=p_c2_input(ismember(p_c2_input,conf_ID_DM));
p_c1_DM=p_c1_input(ismember(p_c2_input,conf_ID_DM));
p_t1_DM=p_t1_input(ismember(p_t1_input,p_c1_DM));
p_t2_DM=p_t2_input(ismember(p_t1_input,p_c1_DM));
p_a1_DM=p_a1_input(ismember(p_a1_input,p_c1_DM));
p_a2_DM=p_a2_input(ismember(p_a1_input,p_c1_DM));


conf_ID_AI=conf_ID_input(conf_label_input==2);
conf_label_AI=2*ones(length(conf_ID_AI),1);
conflist_AI=conflist_input(conf_label_input==2);
p_c2_AI=p_c2_input(ismember(p_c2_input,conf_ID_AI));
p_c1_AI=p_c1_input(ismember(p_c2_input,conf_ID_AI));
p_t1_AI=p_t1_input(ismember(p_t1_input,p_c1_AI));
p_t2_AI=p_t2_input(ismember(p_t1_input,p_c1_AI));
p_a1_AI=p_a1_input(ismember(p_a1_input,p_c1_AI));
p_a2_AI=p_a2_input(ismember(p_a1_input,p_c1_AI));

conf_ID_IR=conf_ID_input(conf_label_input==3);
conf_label_IR=3*ones(length(conf_ID_IR),1);
conflist_IR=conflist_input(conf_label_input==3);
p_c2_IR=p_c2_input(ismember(p_c2_input,conf_ID_IR));
p_c1_IR=p_c1_input(ismember(p_c2_input,conf_ID_IR));
p_t1_IR=p_t1_input(ismember(p_t1_input,p_c1_IR));
p_t2_IR=p_t2_input(ismember(p_t1_input,p_c1_IR));
p_a1_IR=p_a1_input(ismember(p_a1_input,p_c1_IR));
p_a2_IR=p_a2_input(ismember(p_a1_input,p_c1_IR));

% field_included={'DB','DM','AI','IR'};
% 
% expr_conf_ID=[];
% for field=field_included
%     field1=field{1};
%     expr_conf_ID=[expr_conf_ID,'conf_ID_',field1]

conf_ID=[conf_ID_DB;conf_ID_AI;conf_ID_IR];
conflist=[conflist_DB;conflist_AI;conflist_IR];
conf_label=[conf_label_DB;conf_label_AI;conf_label_IR];
p_c2=[p_c2_DB;p_c2_AI;p_c2_IR];
p_c1=[p_c1_DB;p_c1_AI;p_c1_IR];
p_t2=[p_t2_DB;p_t2_AI;p_t2_IR];
p_t1=[p_t1_DB;p_t1_AI;p_t1_IR];
p_a2=[p_a2_DB;p_a2_AI;p_a2_IR];
p_a1=[p_a1_DB;p_a1_AI;p_a1_IR];

termlist_ID=termlist_ID_input(ismember(termlist_ID_input,p_t2));
termlist=termlist_input(ismember(termlist_ID_input,p_t2));
authorlist_ID=authorlist_ID_input(ismember(authorlist_ID_input,p_a2));
authorlist_label=authorlist_label_input(ismember(authorlist_ID_input,p_a2));
authorlist_name=authorlist_name_input(ismember(authorlist_ID_input,p_a2));







for i=1:length(p_a2)
    if(i==1)
        index=1;
    end
    if(isempty(find(authorlist_ID==p_a2(index))))
        p_a2(index)=[];
        p_a1(index)=[];
        index=index-1;
    else
        p_a2(index)=find(authorlist_ID==p_a2(index));
    end
    index=index+1;
end

for i=1:length(p_t2)
    if(i==1)
        index=1;
    end
    if(isempty(find(termlist_ID==p_t2(index))))
        p_t2(index)=[];
        p_t1(index)=[];
        index=index-1;
    else
        p_t2(index)=find(termlist_ID==p_t2(index));
    end
    index=index+1;
end


 for i=1:length(p_c2)
    if(i==1)
        index=1;
    end
    if(isempty(find(conf_ID==p_c2(index))))
        p_c2(index)=[];
        p_c1(index)=[];
        index=index-1;
    else
        p_c2(index)=find(conf_ID==p_c2(index));
    end
    index=index+1;
end
  %remove those in p_a2 or p_t2 or p_c2 which are not contained in authorlist_ID
 %or termlist_ID
 
%p_a2 is the vector whose entries represent the indices for authorlist_name
%or authorlist_label or authorlist_ID



W=sparse([],[],[],max(p_a2),max(p_t2)); %author--term matrix
A0_a_a=sparse([],[],[],max(p_a2),max(p_a2));
A1_t_t=sparse([],[],[],max(p_t2),max(p_t2));
F_t_c=sparse([],[],[],max(p_t2),max(p_c2));
W_a_c=sparse([],[],[],max(p_a2),max(p_c2));

paperlist_full=min(min(p_a1),min(min(p_t1),min(p_c1))):max(max(p_a1),max(max(p_t1),max(p_c1)));
paperlist=paperlist_full(ismember(paperlist_full,[p_t1;p_a1;p_c1]));

for j=1:length(paperlist)
    i=paperlist(j);
    W(p_a2(p_a1==i),p_t2(p_t1==i))=W(p_a2(p_a1==i),p_t2(p_t1==i))+1;
    A0_a_a(p_a2(p_a1==i),p_a2(p_a1==i))=A0_a_a(p_a2(p_a1==i),p_a2(p_a1==i))+1;
    A1_t_t(p_t2(p_t1==i),p_t2(p_t1==i))=A1_t_t(p_t2(p_t1==i),p_t2(p_t1==i))+1;
    F_t_c(p_t2(p_t1==i),p_c2(p_c1==i))=F_t_c(p_t2(p_t1==i),p_c2(p_c1==i))+1;
    W_a_c(p_a2(p_a1==i),p_c2(p_c1==i))=W_a_c(p_a2(p_a1==i),p_c2(p_c1==i))+1;
end
% save('DBLP_labeled.mat');

F_t_a=W';   %F_t_a is the term-author matrix


% 
% for i=1:max(max(p_t1),max(p_c1))
%     
% end



% for i=1:max(max(p_a1),max(p_c1))
%     
% end

%% the following is to remove stopwords 
% clear;
% load('DBLP_labeled_feature.mat');
stopwordslist=textread('stopwords_comma.txt','%s','delimiter',',');
stopwordslist{1}='a';

i_ind=0;
for i=1: length(stopwordslist) 
    if(~isempty(strmatch(stopwordslist{i},termlist,'exact')))
        i_ind=i_ind+1;
        ind_remove(i_ind)=strmatch(stopwordslist{i},termlist,'exact'); 
    end
end

termlist(ind_remove)=[];
termlist_ID(ind_remove)=[];
F_t_a(ind_remove,:)=[];
F_t_c(ind_remove,:)=[];

save('DBLP_labeled_feature_removeDMandSTOP.mat');