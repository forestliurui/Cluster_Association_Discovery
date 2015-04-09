function [TPR, FPR, AUC]=roc_curve(B, S0)
%S0 is the binary groundtruth
    TPR = [];
    FPR = [];

    for th = linspace(-eps,max(abs(B(:)))+eps,1000)
        S = abs(B) > th;
        TPR = [TPR;sum(S(:)&S0(:))/(sum(S0(:))+eps)];
        FPR = [FPR;sum(S(:)&~S0(:))/(sum(~S0(:))+eps)];
    end

    AUC = sum((FPR(1:end-1)-FPR(2:end)).*(TPR(1:end-1)+TPR(2:end))/2);


end
