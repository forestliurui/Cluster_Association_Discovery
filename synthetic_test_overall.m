%test on the synthetic dataset
clear;
iteration=0;
for intensity_noise=0:0.1:0.9
    iteration=iteration+1;
    for j=1:40
        result_4method_temp(j,:)=synthetic_test_orthogonality_function(intensity_noise);
   
    end
    result_4method(iteration,:)=mean(result_4method_temp,1);
end

figure;
hold on;

plot(0:0.1:0.9,result_4method(:,1),'-r*','LineWidth',2);
plot(0:0.1:0.9,result_4method(:,2),'--b+');
plot(0:0.1:0.9,result_4method(:,3),'--ko');
plot(0:0.1:0.9,result_4method(:,4),'--cd');
xlabel('intensity of noise');
ylabel('CAN');

legend('NMTFOC','NMTF\_Chris','DNMTF','RCC');
% legend('NMTFOC','NMTF\_Chris','DNMTF');

hold off;