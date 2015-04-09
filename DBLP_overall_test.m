clear;
iteration=0;
for i=0:0.1:0.9
    iteration=iteration+1;
    for j=1:20
        [result_temp(j,:),MI_conf_temp(j,:),MI_author_temp(j,:),MI_matrix_temp(j,:)]=DBLP_labeled_feature_test_function(i);
    end
    F_result(iteration,:)=mean(result_temp,1);
    MI_conf(iteration,:)=mean(MI_conf_temp);
     MI_author(iteration,:)=mean(MI_author_temp);
     MI_matrix(iteration,:)=mean(MI_matrix_temp);
end
% association
figure;
hold on;

% plot(1:size(data_plot,1),data_plot(:,1),'-b*','LineWidth',2);
%      plot(1:size(data_plot,1),data_plot(:,2),'-rx');
%      plot(1:size(data_plot,1),data_plot(:,3),'-k+');
%    plot(1:size(data_plot,1),data_plot(:,4),'-cs'); %include NMTF_BO
%      plot(1:size(data_plot,1),0.05*ones(size(data_plot,1),1),'g--');

plot(0:0.1:0.9,F_result(:,1),'-r*','LineWidth',2);
plot(0:0.1:0.9,F_result(:,2),'--b+');
plot(0:0.1:0.9,F_result(:,3),'--ko');
plot(0:0.1:0.9,F_result(:,4),'--cd');
xlabel('intensity of noise');
ylabel('CAN');

legend('NMTFOC','NMTF\_Chris','DNMTF','RCC');
% legend('NMTFOC','NMTF\_Chris','DNMTF');

hold off;

% MI conf
figure;
hold on;

plot(0:0.1:0.9,MI_conf(:,1),'-r*','LineWidth',2);
plot(0:0.1:0.9,MI_conf(:,2),'--b+');
plot(0:0.1:0.9,MI_conf(:,3),'--ko');
plot(0:0.1:0.9,MI_conf(:,4),'--cd');
xlabel('intensity of noise');
ylabel('normalized MI');

legend('NMTFOC','NMTF\_Chris','DNMTF','RCC');
% legend('NMTFOC','NMTF\_Chris','DNMTF');

hold off;


% MI author
figure;
hold on;

plot(0:0.1:0.9,MI_author(:,1),'-r*','LineWidth',2);
plot(0:0.1:0.9,MI_author(:,2),'--b+');
plot(0:0.1:0.9,MI_author(:,3),'--ko');
plot(0:0.1:0.9,MI_author(:,4),'--cd');
xlabel('intensity of noise');
ylabel('normalized MI');

legend('NMTFOC','NMTF\_Chris','DNMTF','RCC');
% legend('NMTFOC','NMTF\_Chris','DNMTF');

hold off;

% MI for W matrix
figure;
hold on;

plot(0:0.1:0.9,MI_matrix(:,1),'-r*','LineWidth',2);
plot(0:0.1:0.9,MI_matrix(:,2),'--b+');
plot(0:0.1:0.9,MI_matrix(:,3),'--ko');
plot(0:0.1:0.9,MI_matrix(:,4),'--cd');
xlabel('intensity of noise');
ylabel('MI');

legend('NMTFOC','NMTF\_Chris','DNMTF','RCC');
% legend('NMTFOC','NMTF\_Chris','DNMTF');

hold off;