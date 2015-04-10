%disease gene network
clear;

%eliminate genes that are not in 200 pathways
filename_GO='gene_disease_data\gene_GO\GO_pathway.txt';
ffid_GO=fopen(filename_GO);

tline=fgetl(ffid_GO);


j=1;
while ischar(tline)
    
    GO_pathway{j,:}=strread(tline,'%s','delimiter','\t');
    tline=fgetl(ffid_GO);
    
    j=j+1;
end
fclose(ffid_GO);

list_ind=1;
for i=1:200
    for j=1:length(GO_pathway{i,1})
        if(isempty( strmatch( GO_pathway{i,1}{j},'','exact' ))==1 )
            GO_pathway_genelist{list_ind}=GO_pathway{i,1}{j};
            gene_pathway_relationship(list_ind,i)=1;
            list_ind=list_ind+1;
        else
            break;
        end
    end
end

[max_value_gene,gene_class]=max(gene_pathway_relationship,[],2);

%eliminate genes




filename_pheno_class='gene_disease_data\disease_class.txt';
[disorder_name_raw, gene_symbols, OMIM_ID,chromosome, disorder_class]=textread(filename_pheno_class, '%s%s%d%s%s','delimiter','\t');

% filename_pheno_class='disease_class_tab2.txt';
% [disorder_name_raw, disorder_class, disease_size, disease_degree,diseae_class_degree]=textread(filename_pheno_class, '%s%s%d%d%s','delimiter','\t');


load('gene_disease_data\phenotype_network.mat');
load('gene_disease_data\g_p_network.mat');
load('gene_disease_data\ppi_network.mat');
phenotype_name_raw=phenotype_name;

phenotype_name=lower(phenotype_name_raw);

for i=1:length(disorder_name_raw)
    
    disorder_name{i,1}=lower(disorder_name_raw{i}(2:end-4));
    
    
end
multiple_class=0;
for i=1: length(phenotype_name)
    return_temp=strmatch( phenotype_name{i} ,disorder_name);
    
    if(length(return_temp)>1)
        multiple_class=multiple_class+1;
    end
    
    if(isempty(return_temp))
        substract_pheno_indicator(i)=0;
        substract_pheno_class{i}=[];
    else
        substract_pheno_indicator(i)=1;
        substract_pheno_class{i}=disorder_class{return_temp(1)};
        
    end
end

j=1;
% gene_name_sublist      %indicating the name of genes in the later experiments
%  substract_gene_class     %indicating the pathway each gene belongs to
for i=1: length(gene_name)
    return_temp=strmatch( gene_name{i} ,GO_pathway_genelist);
    
    if(isempty(return_temp))
        substract_gene_indicator(i)=0;
%         substract_gene_class{i}=[];
    else
        substract_gene_indicator(i)=1;
         substract_gene_class_temp(j)=find(gene_pathway_relationship(return_temp(1),:) );
        j=j+1;
    end
end

unique_gene_class_index=unique(sort(substract_gene_class_temp));
for i=1:length(unique_gene_class_index)
        substract_gene_class(substract_gene_class_temp==unique_gene_class_index(i)) =i;

end

for i=1:length(unique(substract_gene_class))
   substract_gene_soft_class{i}=find(substract_gene_class==i); 
   len_gene_soft_class(i)=length(substract_gene_soft_class{i});
   aver_len_gene_soft_class=mean(len_gene_soft_class);
   median_len_gene_soft_class=median(len_gene_soft_class);
end


substract_pheno_indicator=logical(substract_pheno_indicator);
substract_gene_indicator=logical(substract_gene_indicator);
phenotype_network_sub=phenotype_network(substract_pheno_indicator,substract_pheno_indicator);
% g_p_network_sub=g_p_network(:,substract_pheno_indicator);
g_p_network_sub=g_p_network(substract_gene_indicator,substract_pheno_indicator);
ppi_network_sub=ppi_network(substract_gene_indicator,substract_gene_indicator);


j=1;
for i=1:length(substract_pheno_indicator)
    if(substract_pheno_indicator(i)==1)
        substract_pheno_class_sub{j}=substract_pheno_class{i};
        j=j+1;
    end
end

j=1;
pheno_class_list{1,1}='NULL';
for i=1:length(substract_pheno_class_sub)
    if(isempty(strmatch(substract_pheno_class_sub{i}, pheno_class_list)   ) )
        pheno_class_list{j}=substract_pheno_class_sub{i};
        j=j+1;
    end
end

for i=1:length(substract_pheno_class_sub)
    pheno_class_label(i,1)=strmatch(substract_pheno_class_sub{i}, pheno_class_list);  
end

A1=ppi_network_sub;
A2=phenotype_network_sub;
W_a_c=g_p_network_sub;
k1=length(unique_gene_class_index);
k2=21;

iteration_index=0;
% for method_iteration={'NMTFOC','NMTF_BO','DNMTF','RCC'}
for method_iteration={'NMTFOC','NMTF_BO','DNMTF','RCC'}
% for method_iteration={'NMTFOC'}
    iteration_index=iteration_index+1;
    method=method_iteration{1,1};
    switch method
        case 'NMTFOC_ft'
            [ind_H1, ind_H2, S, H1, H2, W]=NMTFOC_ft(A1,A2,F_t_a,F_t_c,k1,k2,W_a_c);
        
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

for i=1:21
    pheno_class_truth4MI{i}=find(pheno_class_label==i);
    pheno_class_predicted4MI{i}=find(ind_H2==(i));
     class_true_label4each_cluster{i}=pheno_class_label(find(ind_H2==(i)));       
end

for i=1:k1
      gene_class_truth4MI{i}=find(substract_gene_class==i);
      gene_class_predicted4MI{i}=find(ind_H1==(i));     
end

% %use mean value of groundtruth cluster size to choose the threshold i
%  for i=min(min(H2)): ( max(max(H2))-min(min(H2)) )/1000 :max(max(H2))
%         H2_multi=(H2>i);
%         if(mean(sum(H2_multi,1)) < aver_len_gene_soft_class  )
%             for j=1:size(H2,2)
%                 gene_multi_class_predicted{j}=find(H2_multi(:,j));
%             end
%             break;
%         end
%  end
 
 
 %use median value of groundtruth cluster size to choose the threshold i
 for i=min(min(H2)): ( max(max(H2))-min(min(H2)) )/1000 :max(max(H2))
        H2_multi=(H2>i);
        if(median(sum(H2_multi,1)) < median_len_gene_soft_class  )
            for j=1:size(H2,2)
                gene_multi_class_predicted{j}=find(H2_multi(:,j));
            end
            break;
        end
 end
 

MI_pheno(iteration_index)=mutual_information_metric(pheno_class_predicted4MI,pheno_class_truth4MI);
MI_gene(iteration_index)=mutual_information_metric(gene_class_predicted4MI,gene_class_truth4MI);
CAM_pheno(iteration_index)=clustering_association_measure(pheno_class_predicted4MI,{},pheno_class_truth4MI,{} );
CAM_gene(iteration_index)=clustering_association_measure(gene_class_predicted4MI,{},gene_class_truth4MI,{} );
CAM_gene_multi(iteration_index)=clustering_association_measure(gene_multi_class_predicted,{},substract_gene_soft_class,{} );

% MI_matrix(iteration_index)=mutual_information_matrix(W_a_c,ind_H1,ind_H2);
fprintf(1,'MI phenotype: %f\n',MI_pheno(iteration_index));
fprintf(1,'MI gene: %f\n',MI_gene(iteration_index));
end
