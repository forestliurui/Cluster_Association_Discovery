function [F_u, G_u, H_u]=RCC(X, WF, WG, k0, k1)
%X is the data mapping matrix 
%WF is the adjacency matrix for graph F
%WG is the adjacency matrix for graph G
%k0 is the number of clusters in graph F
%k1 is the number of clusters in graph G
%X-FHG'-S

[M0, M1]=size(X);
lambda_S=1;
lambda_F=1;
lambda_G=1;

threshold=0.01;
iteration=1;
Iter_threshold=10;


%initialization
F=10*rand(M0, k0);
G=10*rand(M1, k1);
H=10*rand(k0, k1);

%% the beginning of first iteration
    E=X-F*H*(G');
    S_u=(E-(lambda_S/2)*sign(E)).*(abs(E)>(lambda_S/2)); %update S
    S=S_u;
    %some prelimilary work for updating F
    temp_F=[];
    for i=1:M0
        temp_F=[temp_F,(sum((F-ones(M0,1)*F(i,:)).^(2),2)).^(1/2)];
    end
    WF_tilde=WF./(2*temp_F);
    for i=1:M0
        WF_tilde(i,i)=0;
    end
    DF_tilde=diag(sum(WF_tilde,2));
    P=(X-S_u)*G*H';
    Q=H*(G')*G*H';
    
    P_p=(abs(P)+P)./2;
    P_n=(abs(P)-P)./2;
    Q_p=(abs(Q)+Q)./2;
    Q_n=(abs(Q)-Q)./2;
    A=(P_n+F*Q+lambda_F*DF_tilde*F)./F;
    C=(P_p+lambda_F*WF_tilde*F).*F;
    F_u=Algo_interior(A, C, F); %update F
    
    H_numerator=(F_u')*(X-S_u)*G;
    H_u=H.*((((abs(H_numerator)+H_numerator)./2)./((F_u')*F*H*(G')*G+(abs(H_numerator)-H_numerator)./2)).^(1/2)); %update H
    
    %some preliminary work for updating G
    temp_G=[];
    for i=1:M1
        temp_G=[temp_G,(sum((G-ones(M1,1)*G(i,:)).^(2),2)).^(1/2)];
    end
    WG_tilde=WG./(2*temp_G);
    for i=1:M1
        WG_tilde(i,i)=0;
    end
    DG_tilde=diag(sum(WG_tilde,2));
    P=((X-S_u)')*F_u*H_u;
    Q=(H_u')*(F_u')*F_u*H_u;
    
    P_p=(abs(P)+P)./2;
    P_n=(abs(P)-P)./2;
    Q_p=(abs(Q)+Q)./2;
    Q_n=(abs(Q)-Q)./2;
    A=(P_n+G*Q+lambda_G*DG_tilde*G)./G;
    C=(P_p+lambda_G*WG_tilde*G).*G;
    G_u=Algo_interior(A, C, G); %update G
    
    %compute residue
    S_res=norm(S_u-S);
    F_res=norm(F_u-F);   %F corresponds to H0
   G_res=norm(G_u-G);     %G corresponds to H1
%% the end of first iteration

while(((S_res>threshold)||(F_res>threshold)||(G_res>threshold))&&(iteration<Iter_threshold))
% while(iteration<Iter_threshold)
    iteration=iteration+1;
    S=S_u;
    F=F_u;
    G=G_u;
    
    E=X-F*H*(G');
    S_u=(E-(lambda_S/2)*sign(E)).*(abs(E)>(lambda_S/2)); %update S
    
    %some prelimilary work for updating F
    temp_F=[];
    for i=1:M0
        temp_F=[temp_F,(sum((F-ones(M0,1)*F(i,:)).^(2),2)).^(1/2)];
    end
    WF_tilde=WF./(2*temp_F);
    for i=1:M0
        WF_tilde(i,i)=0;
    end
    DF_tilde=diag(sum(WF_tilde,2));
    P=(X-S_u)*G*H';
    Q=H*(G')*G*H';
    
    P_p=(abs(P)+P)./2;
    P_n=(abs(P)-P)./2;
    Q_p=(abs(Q)+Q)./2;
    Q_n=(abs(Q)-Q)./2;
    A=(P_n+F*Q+lambda_F*DF_tilde*F)./F;
    C=(P_p+lambda_F*WF_tilde*F).*F;
    F_u=Algo_interior(A, C, F); %update F
    
    H_numerator=(F_u')*(X-S_u)*G;
    H_u=H.*((((abs(H_numerator)+H_numerator)./2)./((F_u')*F*H*(G')*G+(abs(H_numerator)-H_numerator)./2)).^(1/2)); %update H
    
    %some preliminary work for updating G
    temp_G=[];
    for i=1:M1
        temp_G=[temp_G,(sum((G-ones(M1,1)*G(i,:)).^(2),2)).^(1/2)];
    end
     WG_tilde=WG./(2*temp_G);
    for i=1:M1
        WG_tilde(i,i)=0;
    end
    DG_tilde=diag(sum(WG_tilde,2));
    P=((X-S_u)')*F_u*H_u;
    Q=(H_u')*(F_u')*F_u*H_u;
    
    P_p=(abs(P)+P)./2;
    P_n=(abs(P)-P)./2;
    Q_p=(abs(Q)+Q)./2;
    Q_n=(abs(Q)-Q)./2;
    A=(P_n+G*Q+lambda_G*DG_tilde*G)./G;
    C=(P_p+lambda_G*WG_tilde*G).*G;
    G_u=Algo_interior(A, C, G); %update G
    
    %compute residue
    S_res=norm(S_u-S);
    F_res=norm(F_u-F);
    G_res=norm(G_u-G);
end

% pause;
    
    
    