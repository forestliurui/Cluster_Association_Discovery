function [ind_H0, ind_H1, S, H0, H1, W]=NMTFOC_ft(varargin)
%NMTFOC with association matrix learnt from the feature matrix
%the formular is as follows
%norm(W_01-H0*S01*H1')+alpha*Tr(H0'*Theta_0*H0)+beta*Tr(H1'*Theta_1*H1)+eta1*norm(S01,1)+gamma*norm(F1-F0*W01)+eta0*norm(W01,1)

A0=varargin{1};
A1=varargin{2};
F0=varargin{3};
F1=varargin{4};
k0=varargin{5};
k1=varargin{6};




alpha=1;  
beta=1;
eta1=1;
threshold=0.01;
iteration=1;
eta0=20;
gamma=1;

[M0,N0]=size(A0);
[M1,N1]=size(A1);

%initialization of W
if(nargin==7)
    W=varargin{7};
else
    W=10*rand(N0,M1);
end

Theta0=diag(sum(A0,2))-A0;
Theta1=diag(sum(A1,2))-A1;

H0=10*rand(M0, k0);
H1=10*rand(M1, k1);
S=10*rand(k0,k1);
% S(1,1)=30;
% S(2,2)=30;
% W_p=(abs(W)+W)./2;
% W_n=(abs(W)-W)./2;
Theta0_n=(abs(Theta0)-Theta0)./2;
Theta0_p=(abs(Theta0)+Theta0)./2;
Theta1_n=(abs(Theta1)-Theta1)./2;
Theta1_p=(abs(Theta1)+Theta1)./2;
F01_p=(abs((F0')*F1)+(F0')*F1)./2;
F01_n=(abs((F0')*F1)-(F0')*F1)./2;
F00_p=(abs((F0')*F0)+(F0')*F0)./2;
F00_n=(abs((F0')*F0)-(F0')*F0)./2;
W_n=10*zeros(N0,M1);


W_u=W.*((H0*S*(H1')+(gamma*F01_p)+0.5*(F00_n)*W)./(W+gamma*F01_n+gamma*F00_p*W+0.5*eta0)).^(1/2);
S_u=S.*((H0'*W*H1)./((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )).^(1/2);
H0_u=H0.*((W*H1*S'+alpha*Theta0_n*H0)./((H0*H0'*(W*H1*S'+alpha*Theta0_n*H0)))).^(1/2);
H1_u=H1.*((W'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W'*H0*S+beta*Theta1_n*H1) )).^(1/2);

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
    W=W_u;
%     S_n= (((H0'*W*H1)-0.5*eta1)./((H0')*H0*S*(H1')*H1)).*S;
%     H0_n= ((W*H1*(S'))./(alpha*(Theta0*H0)+(H0*S*(H1')*H1*(S')  ) )).*H0;
%     
%     H0_n= H0_n./(sum(H0_n,2)*ones(1,k0));
%     
%     H1_n= (((W')*H0*S)./(beta*(Theta1*H1)+(H1*(S')*(H0')*H0*S  ) )).*H1;
%     
%     H1_n= H1_n./(sum(H1_n,2)*ones(1,k1));

    W_u=W.*((H0*S*(H1')+(gamma*F01_p)+0.5*(F00_n)*W)./(W+gamma*F01_n+gamma*F00_p*W+0.5*eta0)).^(1/2);
    S_u=S.*((H0'*W*H1)./((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )).^(1/2);
    
    H0_u=H0.*((W*H1*S'+alpha*Theta0_n*H0)./((H0*H0'*(W*H1*S'+alpha*Theta0_n*H0)))).^(1/2);    
    H0_u= H0_u./(sum(H0_u,2)*ones(1,k0));
    
    H1_u=H1.*((W'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W'*H0*S+beta*Theta1_n*H1) )).^(1/2);      
    H1_u= H1_u./(sum(H1_u,2)*ones(1,k1));
    
    S_res=norm(S_u-S);
    H0_res=norm(H0_u-H0);
    
    H1_res=norm(H1_u-H1);
end

[val_H0, ind_H0]=max(H0,[],2);
[val_H1, ind_H1]=max(H1,[],2);

f_res=norm(W-(H0)*S*(H1'))+alpha*trace((H0')*Theta0*H0)+beta*trace((H1')*Theta1*H1)+eta1*norm(S(:),1)
% pause;

