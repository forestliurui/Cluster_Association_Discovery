clear;
% 
% [SNP_chr	SNP_start1	SNP_end1	SNP_n	SNP_start2	SNP_end2	SNP_n2	SNP_pos]=textread('07SNP_information_1017_ordered_expe.txt','%d%d%d%d%s%s%s%d');
% 
% [Gene_ID	Gene_ORF	Gene_geneSym	Gene_chr Gene_end Gene_start]=textread('06Gene_information_4474_ordered_expe.txt','%d%s%s%d%d%d');
% 
% 
% iteration=1;
% for i=1:length(SNP_start1)
%     temp=find(Gene_chr==SNP_chr(i));
%     for t=1:length(temp)
%         j=temp(t);
%         if(( (  SNP_start1(i)>Gene_start(j))&& (  SNP_start1(i)<Gene_end(j) ))||((SNP_end1(i)>Gene_start(j))&&(SNP_end1(i)<Gene_end(j)))||((SNP_start1(i)<Gene_end(j))&&(SNP_end1(i)>Gene_end(j))))
%             SNP_mapping(iteration)=i;
%             Gene_mapping(iteration)=j;
%             iteration=iteration+1;
%         end
%     end
% end
% save('SNP_Gene_yeat_mapping.mat');

load('SNP_Gene_yeat_mapping.mat');

%load('data_output\results_yeast_NMTFMG_2014_9_2_11_51.mat'); %to process
%NMTFMG data
% load('data_output\results_yeast_NMTFOC_2014_9_3_10_43.mat');  %to process
% NMTFOC data
load('data\results_yeast_NMTFOC_2014_10_6_15_39');

[Gene_ID	Gene_ORF	Gene_geneSym	Gene_chr_te Gene_end_te Gene_start_te]=textread('06Gene_information_4474_ordered_expe1.txt','%d%s%s%s%s%s');

filename=['results_yeast_gene_cluster_NMTFOC.txt'];
ffile=fopen(filename,'wt');
S_temp=S;   
for i=1:10
[rows, cols]=find(max(max(S_temp))==S_temp);


[SNP_map_v,SNP_map_ind_temp]=ismember((find(ind_H1==rows)),SNP_mapping);
SNP_map_ind=SNP_map_ind_temp(find(SNP_map_ind_temp~=0));

index_union_SNPmapping_Gene=union(Gene_mapping(SNP_map_ind),find(ind_H2==cols));
index_SNPmapping=Gene_mapping(SNP_map_ind);
index_Gene=find(ind_H2==cols);



fprintf(ffile,'%s\t',Gene_geneSym{index_SNPmapping});
fprintf(ffile,'\n');

fprintf(ffile,'%s\t',Gene_geneSym{index_Gene});
fprintf(ffile,'\n');

fprintf(ffile,'%s\t',Gene_geneSym{index_union_SNPmapping_Gene});
fprintf(ffile,'\n');
%fprintf(ffile,'the majority of authors from author cluster 
%d are from %d\n\n\n',i,voting_strategy(authorlist_label(find(ind_H1==rows))+1)-1);

S_temp(rows,cols)=0;
end
 
fclose(ffile);



filename=['results_yeast_gene_cluster_NMTFOC_onlyGene.txt'];
ffile=fopen(filename,'wt');
 
for i=1:size(H2)

index_Gene=find(ind_H2==i);

fprintf(ffile,'%s\t',Gene_geneSym{index_Gene});
fprintf(ffile,'\n');

%fprintf(ffile,'the majority of authors from author cluster 
%d are from %d\n\n\n',i,voting_strategy(authorlist_label(find(ind_H1==rows))+1)-1);

end
 
fclose(ffile);


filename=['results_yeast_gene_cluster_NMTFOC_onlySNP.txt'];
ffile=fopen(filename,'wt');
 
for i=1:size(H1)

[SNP_map_v,SNP_map_ind_temp]=ismember((find(ind_H1==i)),SNP_mapping);
SNP_map_ind=SNP_map_ind_temp(find(SNP_map_ind_temp~=0));
index_SNPmapping=Gene_mapping(SNP_map_ind);


fprintf(ffile,'%s\t',Gene_geneSym{index_SNPmapping});
fprintf(ffile,'\n');

%fprintf(ffile,'the majority of authors from author cluster 
%d are from %d\n\n\n',i,voting_strategy(authorlist_label(find(ind_H1==rows))+1)-1);

end
 
fclose(ffile);

