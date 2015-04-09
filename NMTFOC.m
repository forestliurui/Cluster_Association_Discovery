function [ind_H0, ind_H1, S, H0, H1]=NMTFOC(W, A0, A1, k0, k1)

alpha=1; 
beta=1;
eta1=1;
threshold=0.01;
iteration=1;

[M0,N0]=size(A0);
[M1,N1]=size(A1);

Theta0=diag(sum(A0,2))-A0;
Theta1=diag(sum(A1,2))-A1;


% seed = 2345;               % set the seed for random number generator   
% rng(seed);
H0=10*rand(M0, k0);
H1=10*rand(M1, k1);
S=10*rand(k0,k1);
% S(1,1)=30;
% S(2,2)=30;
W_p=(abs(W)+W)./2;
W_n=(abs(W)-W)./2;
Theta0_n=(abs(Theta0)-Theta0)./2;
Theta0_p=(abs(Theta0)+Theta0)./2;
Theta1_n=(abs(Theta1)-Theta1)./2;
Theta1_p=(abs(Theta1)+Theta1)./2;

S_u=S.*((H0'*W_p*H1)./((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )).^(1/2);
H0_u=H0.*((W_p*H1*S'+alpha*Theta0_n*H0)./((H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0)))).^(1/2);
H1_u=H1.*((W_p'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1) )).^(1/2);

%         S_u=S.* ( gdivide_revised((H0'*W_p*H1) ,(H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1   )).^(1/2);
%         H0_u=H0.*(gdivide_revised(W_p*H1*S'+alpha*Theta0_n*H0 ,  (H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0))   )).^(1/2);
%          H0_u= H0_u./(sum(H0_u,2)*ones(1,k0));
%          H1_u=H1.*(gdivide_revised((W_p'*H0*S+beta*Theta1_n*H1),H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1)  )   ).^(1/2)  ;
%          H1_u= H1_u./(sum(H1_u,2)*ones(1,k1));



if(sum(sum(isnan(S_u)))>0)
    ;
end
if(sum(sum(isnan(H0_u)))>0)
    ;
end
if(sum(sum(isnan(H1_u)))>0)
    ;
end

% S_n= (((H0'*W*H1)-0.5*eta1)./(H0'*H0*S*(H1')*H1)).*S;
% H0_n= (W*H1*(S'))./(alpha*(Theta0*H0)+(H0*S*(H1')*H1*S'  ) ).*H0;
% H1_n= ((W')*H0*S)./(beta*(Theta1*H1)+(H1*(S')*(H0')*H0*S  ) ).*H1;

S_res=norm(S_u-S);
H0_res=norm(H0_u-H0);
H1_res=norm(H1_u-H1);

while((S_res>threshold)||(H0_res>threshold)||(H1_res>threshold))
    iteration=iteration+1;
    S=S_u;
    H0=H0_u;
    H1=H1_u;
%     S_n= (((H0'*W*H1)-0.5*eta1)./((H0')*H0*S*(H1')*H1)).*S;
%     H0_n= ((W*H1*(S'))./(alpha*(Theta0*H0)+(H0*S*(H1')*H1*(S')  ) )).*H0;
%     
%     H0_n= H0_n./(sum(H0_n,2)*ones(1,k0));
%     
%     H1_n= (((W')*H0*S)./(beta*(Theta1*H1)+(H1*(S')*(H0')*H0*S  ) )).*H1;
%     
%     H1_n= H1_n./(sum(H1_n,2)*ones(1,k1));



    S_u=S.*(  (H0'*W_p*H1)    ./   ((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )     ).^(1/2);

    
     
   H0_u=H0.*((W_p*H1*S'+alpha*Theta0_n*H0)./(  (H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0))   )  ).^(1/2);
    
    H0_u= H0_u./(sum(H0_u,2)*ones(1,k0));
    
    
    
    
    
    H1_u=H1.*((W_p'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1) )).^(1/2);
    
    H1_u= H1_u./(sum(H1_u,2)*ones(1,k1));
    
    
%         S_u=S.* ( gdivide_revised((H0'*W_p*H1) ,(H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1   )).^(1/2);
%         H0_u=H0.*(gdivide_revised(W_p*H1*S'+alpha*Theta0_n*H0 ,  (H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0))   )).^(1/2);
%          H0_u= H0_u./(sum(H0_u,2)*ones(1,k0));
%          H1_u=H1.*(gdivide_revised((W_p'*H0*S+beta*Theta1_n*H1),H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1)  )   ).^(1/2)  ;
%          H1_u= H1_u./(sum(H1_u,2)*ones(1,k1));
        
    if(sum(sum(isnan(S_u)))>0)
    ;
    end
    if(sum(sum(isnan(H0_u)))>0)
    ;
    end
    if(sum(sum(isnan(H1_u)))>0)
    ;
    end
    
    S_res=norm(S_u-S);
    H0_res=norm(H0_u-H0);
    
    H1_res=norm(H1_u-H1);
end

[val_H0, ind_H0]=max(H0,[],2);
[val_H1, ind_H1]=max(H1,[],2);

f_res=norm(W-(H0)*S*(H1'))+alpha*trace((H0')*Theta0*H0)+beta*trace((H1')*Theta1*H1)+eta1*norm(S(:),1)
% pause;

