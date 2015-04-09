%to compute the calibrated p value for Gene Ontology Analysis on Yeast
clear;
% input_file_directory='k1_100k2_300hardclus_nothreshold';
input_file_directory='yeast_hard_clus_k100300_PPI';
iteration=0;
% for method_iteration={'3HNMTFOCSpec_SNP','NMTF_BO_SNP','DNMTF_SNP','RCC_SNP'}
 for method_iteration={'NMTFOC_gene','NMTF_BO_gene','DNMTF_gene','RCC_gene'}
%      for method_iteration={'NMTFOC_SNP','NMTF_BO_SNP','DNMTF_SNP','RCC_SNP'}
%        for method_iteration={'NMTFOC_genenaturescience','NMTF_BO_genenaturescience'}
%      for method_iteration={'NMTFOC_SNPnaturescience'}
    iteration=iteration+1;
    method=method_iteration{1,1};
    
%    input_file_test_null=['k1_100k2_100_SoftClustering\results_yeast_gene_cluster_',method,'.txt'];
%     input_file_test_null=['k1_100k2_300_HardClustering\results_yeast_gene_cluster_',method,'.txt'];
      input_file_test_null=[input_file_directory,'\results_yeast_gene_cluster_',method,'.txt'];
    test_null_string=textread(input_file_test_null,'%s','delimiter','\n','bufsize',40000);
    
    input_file_GO=[input_file_directory,'\gene_ontology_yeast.txt'];
    [GO_name, GO_description]=textread(input_file_GO, '%s%s','delimiter','\t');
    
%     fid_GO = fopen(input_file_GO);  %open file
%     data_GO = fread(fid_GO, '*char')'; %read all contents into data as a char array (don't forget the `'` to make it a row rather than a column).
%     fclose(fid_GO);
%     entries = regexp(data_GO, ',', 'split');
    
    if(iteration==1)
        input_file_test_association=[input_file_directory,'\results_yeast_gene_cluster_NMTFOC_assotop100.txt'];
    elseif(iteration==2)
        input_file_test_association=[input_file_directory,'\results_yeast_gene_cluster_NMTF_BO_assotop100.txt'];
    elseif(iteration==3)
        input_file_test_association=[input_file_directory,'\results_yeast_gene_cluster_DNMTF_assotop100.txt'];
    else
        input_file_test_association=[input_file_directory,'\results_yeast_gene_cluster_RCC_assotop100.txt'];
    end
    test_association_string=textread(input_file_test_association,'%s','delimiter','\n','bufsize',40000);
    
    j_gene=1;
    j_SNP=1;
    for i=1:length(test_association_string)
        if(mod(i,3)==2)
            test_association_string_gene{j_gene}=test_association_string{i};
            j_gene=j_gene+1;
        elseif(mod(i,3)==1)
            test_association_string_SNP{j_SNP}=test_association_string{i};
            j_SNP=j_SNP+1;
        end
        
    end
    
    
    
%     input_file_RawPValue=['k1_100k2_300_HardClustering\',method,'\S02RawPValue.log'];
    input_file_RawPValue=[input_file_directory,'\',method,'\S02RawPValue.log'];
%         input_file_RawPValue=['yeast_hard_clus_k100300_PPI\',method,'\S02RawPValue.log'];
    %read the Raw PValue
    [clusterID, GO_ID, RawPValue]=textread(input_file_RawPValue,'%d:%d:%f','headerlines',18);  %clusterDI and GO_ID start from 0
    
%     output_file_calibratedPValue=['k1_100k2_300_HardClustering\',method,'\S06CalibratedPValue.log'];
     output_file_calibratedPValue=[input_file_directory,'\',method,'\S06CalibratedPValue.log'];
    %in our current setting, the clusterID is from 0 to 29

for index_permutation=1:length(RawPValue)
    if(isempty(test_null_string{index_permutation})==0)
    index_permutation_file=index_permutation-1;
%     input_file_permutation=['k1_100k2_300_HardClustering\',method,'\S05With100PermNullDist4OneSet',num2str(index_permutation_file),'.log'];
    input_file_permutation=[input_file_directory,'\',method,'\S05With100PermNullDist4OneSet',num2str(index_permutation_file),'.log'];
    [permutationID, temp_perm, PValue_perm]=textread(input_file_permutation,'%d:%d:%f','headerlines',18);
  
    num_permutation=length(PValue_perm);
    
    PValue_perm_Raw=[RawPValue(index_permutation);PValue_perm];
            
    [PValue_perm_sorted,PValue_perm_sorted_index]=sort(PValue_perm_Raw,'ascend');
    
    position_RawPValue=find(PValue_perm_sorted_index==1);
    
    calibrated_PValue(index_permutation)=position_RawPValue/num_permutation;      
    else
        calibrated_PValue(index_permutation)=1;
    end
end

    for i=1:length(test_association_string_gene)
        for j=1:length(test_null_string)
            if(strmatch(test_null_string{j},test_association_string_gene{i},'exact')==1)
                break;
            end
        end
        test_association_string_gene_PValue(i,iteration)=calibrated_PValue(j);
        test_association_string_gene_GO_indice(i,iteration)=GO_ID(j);
    end
    
     for i=1:length(test_association_string_SNP)
        for j=1:length(test_null_string)
            if(strmatch(test_null_string{j},test_association_string_SNP{i},'exact')==1)
                break;
            end
        end
        test_association_string_SNP_PValue(i,iteration)=calibrated_PValue(j);
        test_association_string_SNP_GO_indice(i,iteration)=GO_ID(j);
    end

ffile=fopen(output_file_calibratedPValue,'wt');

fprintf(ffile,'%f\n',calibrated_PValue);
fclose(ffile);


% output for GO association result 

output_file_GO_indices_association=[input_file_directory,'\GO_indices_association_',method,'.txt'];
ffile_output_file_GO_indice=fopen(output_file_GO_indices_association,'wt');
% fprintf(ffile_output_file_GO_indice,'%d\n',test_association_string_SNP_GO_indice);
fprintf(ffile_output_file_GO_indice,'%d\n',test_association_string_gene_GO_indice);
fclose(ffile_output_file_GO_indice);

%output for GO association result



% data_plot_temp=[];
% for i=0:9
%     data_plot_temp=[data_plot_temp;(calibrated_PValue(3*i+1:3*(i+1)))];
% end    
% data_plot(:,iteration)=data_plot_temp;

data_plot_o(:,iteration)=calibrated_PValue;

%     b=bar(data_plot);
%     ch=get(b,'children');
%     set(gca,'XTickLabel',{'Pair 1','Pair 2','Pair 3','Pair 4','Pair 5','Pair 6','Pair 7','Pair 8','Pair 9','Pair 10'});
%     legend('SNP','gene','combined');
%     title(method);
end
% data_plot_o=test_association_string_gene_PValue;
data_plot_o=test_association_string_SNP_PValue;

%%
[data_plot,data_plot_ind]=sort(data_plot_o,1);
% [GO_name, GO_description]

output_GO=[input_file_directory,'\GO_yeast.txt'];
ffile_output_GO=fopen(output_GO,'wt');

fprintf(ffile_output_GO,'%s\n',  GO_name{GO_ID(data_plot_ind(:,1))+1});


% data_plot=data_plot_o;
figure;
hold on;
     plot(1:size(data_plot,1),data_plot(:,1),'-r*','LineWidth',2);
     plot(1:size(data_plot,1),data_plot(:,2),'--b+');
     plot(1:size(data_plot,1),data_plot(:,3),'--ko');
   plot(1:size(data_plot,1),data_plot(:,4),'--cd'); %include NMTF_BO
     plot(1:size(data_plot,1),0.05*ones(size(data_plot,1),1),'g--');
%    legend('NMTFOC\_SNP','DNMTF\_SNP','RCC\_SNP','NMTF\_Chris\_SNP');
    legend('NMTFOC\_gene','NMTF\_Chris\_gene','DNMTF\_gene','RCC\_gene');
%      legend('NMTFOC\_gene','NMTF\_Chris\_gene','DNMTF\_gene','RCC\_gene');
     ylabel('Calibrated P Value');
%      xlabel('set ID');
      xlabel('top set');
%      axis([0 300 0 0.2]);
     axis([0 100 0 0.2]);
hold off;



% 
% %%
% figure;
% hold on;
% %     plot(1:size(data_plot,1),data_plot(:,1,1),'-b*');
% %     plot(1:size(data_plot,1),data_plot(:,2,1),'-bo');
% %     plot(1:size(data_plot,1),data_plot(:,3,1),'-b+');
% %         
% %     plot(1:size(data_plot,1),data_plot(:,1,2),'-k*');
% %     plot(1:size(data_plot,1),data_plot(:,2,2),'-ko');
% %     plot(1:size(data_plot,1),data_plot(:,3,2),'-k+');
% %     
% %     plot(1:size(data_plot,1),data_plot(:,1,3),'-r*');
% %     plot(1:size(data_plot,1),data_plot(:,2,3),'-ro');
% %     plot(1:size(data_plot,1),data_plot(:,3,3),'-r+');
%     
%         plot(1:size(data_plot,1),data_plot(:,1),'-b*');
%      plot(1:size(data_plot,1),data_plot(:,2),'-bo');
%      plot(1:size(data_plot,1),data_plot(:,2),'-b+');
% 
%  hold off;
% 
% %  legend('SNP-NMTFMG','gene-NMTFMG','combined-NMTFMG','SNP-DNMTF','gene-DNMTF','combined-DNMTF','SNP-RCC','gene-RCC','combined-RCC');
%   legend('SNP-NMTFMG','SNP-DNMTF','SNP-RCC');
% %  title(method);
%  xlabel('set');
%  ylabel('p value');

