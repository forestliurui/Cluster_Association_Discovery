function [ind_H,  S, H, W]=NMTFMG_p(varargin)
%A_cell, F_cell,c
%W_cell contains the W's
%A_cell contains the A's
%F_cell contains the F's
%c is the vector who contains the number of clusters for each domain

A_cell=varargin{1};
F_cell=varargin{2};
c=varargin{3};

    
    

display('The beginning of NMTFMG\n');
N=size(A_cell,1);
m=size(F_cell{1,1},1);

for i=1:size(A_cell,2)
    n(i)=size(A_cell{i,i},2);
end


%random initialization of S, H
for i=1:N
    H_cell{i,i}=10*rand(n(i),c(i));
    for j=1:N
        if (i~=j)
            S_cell{j,i}=10*rand(c(j),c(i));
            W_cell{j,i}=10*rand(n(j),n(i));
        end 

    end
end

if(nargin==4)
    W_cell=varargin{4};
end

W=[]; %construct the matirx W
S=[];
I=[];
for i=1:size(W_cell,2)
    W_ctr=[];
    S_ctr=[];
    I_ctr=[];
    for j=1:size(W_cell,1)
        if(i==j)
            W_ctr=[W_ctr;zeros(n(i))];
            S_ctr=[S_ctr;zeros(c(i))];
        else
            W_ctr=[W_ctr;W_cell{j,i}];
            S_ctr=[S_ctr;S_cell{j,i}];
        end
        I_ctr=[I_ctr;eye(m)];
    end
    W=[W,W_ctr];
    S=[S,S_ctr];
    I=[I,I_ctr];
end
W=sparse(W);
S=sparse(S);

W0=W;


alpha=10000*ones(1,N);
% alpha(1)=200000;
for i=1:N
    if(i==1)
        A=alpha(i)*A_cell{1,1};
        F=F_cell{1,1};
        H=H_cell{1,1};
    else
        A=blkdiag(A, alpha(i)*A_cell{i,i});
        F=blkdiag(F,F_cell{i,i});
        H=blkdiag(H, H_cell{i,i}); 
    end
end
A=sparse(A);
F=sparse(F);
H=sparse(H);
Theta=diag(sum(A,2))-A;


% alpha=1; 
% beta=1;
eta1=10;
eta2=0;
eta3=10;
eta4=100;
threshold=0.0001;
iteration=1;
pert=0.0001;
eta0=1000;
% [M0,N0]=size(A0);
% [M1,N1]=size(A1);

% Theta0=diag(sum(A0,2))-A0;
% Theta1=diag(sum(A1,2))-A1;

% H0=10*rand(M0, k0);
% H1=10*rand(M1, k1);
% S=10*rand(k0,k1);
% S(1,1)=30;
% S(2,2)=30;

Theta_n=(abs(Theta)-Theta)./2;
Theta_p=(abs(Theta)+Theta)./2;
% Theta1_n=(abs(Theta1)-Theta1)./2;
% Theta1_p=(abs(Theta1)+Theta1)./2;

W_p=(abs(W)+W)./2;
W_n=(abs(W)-W)./2;

% W_u=(H*S*(H')+eta3*(F')*I*F-eta3*(F')*F)/(1+eta0);  %without nonnegative
% constraint on W

F_double_p=(abs((F')*F)+(F')*F)./2;
F_double_n=(abs((F')*F)-(F')*F)./2;
F_I_p=(abs((F')*I*F)+(F')*I*F)./2;
F_I_n=(abs((F')*I*F)-(F')*I*F)./2;

%W_u=W.*((2*(H*S*(H'))+eta3*(F_double_n*W)+2*eta3*(F_I_p))./(2*W+eta0+2*eta3*(F_double_p*W)+2*eta3*F_I_n)).^(1/2); %impose nonnegative constraint on W

W_u=W.*((2*(H*S*(H'))+eta3*(F_double_n*W)+2*eta3*(F_I_p)+2*eta4*W0)./(2*W+eta0+2*eta3*(F_double_p*W)+2*eta3*F_I_n+2*eta4*W)).^(1/2);


S_u=S.*((((H')*W_p*H)./((H')*W_n*H+((H')*H*S*(H')*H)+0.5*eta1+eta2*S)).^(1/2));
H_u=H.*(((2*W_p*H*S+Theta_n*H)./(H*(H')*(4*W_p*H*S+2*Theta_n*H)+pert*ones(size(H)))).^(1/4));
 H_u= H_u./(sum(H_u,2)*ones(1,size(H,2)));

 
%  H_u=H.*((Theta_n*H+2*eta3*((F')*I*F*H*S))./(H*(H')*(2*Theta_n*H+4*eta3*((F')*I*F*H*S)))).^(1/4);
%  S_u=S.*((2*eta3*((H')*(F')*I*F*H))./(eta1+2*eta2*((H')*(F')*F*H*S*(H')*H))).^(1/2);
%   H_u= H_u./(sum(H_u,2)*ones(1,size(H,2)));

% S_u=S.*((H0'*W_p*H1)./((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )).^(1/2);
% H0_u=H0.*((W_p*H1*S'+alpha*Theta0_n*H0)./((H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0)))).^(1/2);
% H1_u=H1.*((W_p'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1) )).^(1/2);

% S_n= (((H0'*W*H1)-0.5*eta1)./(H0'*H0*S*(H1')*H1)).*S;
% H0_n= (W*H1*(S'))./(alpha*(Theta0*H0)+(H0*S*(H1')*H1*S'  ) ).*H0;
% H1_n= ((W')*H0*S)./(beta*(Theta1*H1)+(H1*(S')*(H0')*H0*S  ) ).*H1;

W_res=norm(full(W_u-W));

S_res=norm(full(S_u-S));
H_res=norm(full(H_u-H));

Iter_threshold=1000;

while(((S_res>threshold)||(H_res>threshold)||(W_res>threshold))&&(iteration<Iter_threshold))
% while((S_res>threshold)||(H_res>threshold)||(W_res>threshold))
    iteration=iteration+1;
    display(['NMTFMG: Iteration ', int2str(iteration),'\n']);
    S=S_u;
    H=H_u;
    W=W_u;
%     S_n= (((H0'*W*H1)-0.5*eta1)./((H0')*H0*S*(H1')*H1)).*S;
%     H0_n= ((W*H1*(S'))./(alpha*(Theta0*H0)+(H0*S*(H1')*H1*(S')  ) )).*H0;
%     
%     H0_n= H0_n./(sum(H0_n,2)*ones(1,k0));
%     
%     H1_n= (((W')*H0*S)./(beta*(Theta1*H1)+(H1*(S')*(H0')*H0*S  ) )).*H1;
%     
%     H1_n= H1_n./(sum(H1_n,2)*ones(1,k1));
% 
%     S_u=S.*((H0'*W_p*H1)./((H0'*W_n*H1)+(H0'*H0*S*(H1')*H1)+0.5*eta1  )).^(1/2);
%     
% %    S_u=S_u/sum(sum(S_u));
%      
%    H0_u=H0.*((W_p*H1*S'+alpha*Theta0_n*H0)./((H0*H0'*(W_p*H1*S'+alpha*Theta0_n*H0)))).^(1/2);
%     
%     H0_u= H0_u./(sum(H0_u,2)*ones(1,k0));
%     
%     H1_u=H1.*((W_p'*H0*S+beta*Theta1_n*H1)./( H1*H1'*(W_p'*H0*S+beta*Theta1_n*H1) )).^(1/2);
%     
%     H1_u= H1_u./(sum(H1_u,2)*ones(1,k1));

    W_p=(abs(W)+W)./2;
    W_n=(abs(W)-W)./2;

%     W_u=(H*S*(H')+eta3*(F')*I*F-eta3*(F')*F)/(1+eta4);
%W_u=W.*((2*(H*S*(H'))+eta3*(F_double_n*W)+2*eta3*(F_I_p))./(2*W+eta0+2*eta3*(F_double_p*W)+2*eta3*F_I_n)).^(1/2);
 
W_u=W.*((2*(H*S*(H'))+eta3*(F_double_n*W)+2*eta3*(F_I_p)+2*eta4*W0)./(2*W+eta0+2*eta3*(F_double_p*W)+2*eta3*F_I_n+2*eta4*W)).^(1/2);


    S_u=S.*((((H')*W_p*H)./((H')*W_n*H+((H')*H*S*(H')*H)+0.5*eta1+eta2*S)).^(1/2));
    H_u=H.*(((2*W_p*H*S+Theta_n*H)./(H*(H')*(4*W_p*H*S+2*Theta_n*H)+pert*ones(size(H)))).^(1/4));
 H_u= H_u./(sum(H_u,2)*ones(1,size(H,2)));

%  S_u=S.*((((H')*W_p*H)./((H')*W_n*H+((H')*H*S*(H')*H)+0.5*eta1+eta2*S)).^(1/2));
%     H_u=H.*(((2*W_p*H*S+Theta_n*H)./(H*(H')*(4*W_p*H*S+2*Theta_n*H)+pert*ones(size(H)))).^(1/4));
%  H_u= H_u./(sum(H_u,2)*ones(1,size(H,2)));


 W_res=norm(full(W_u-W));

S_res=norm(full(S_u-S));
H_res=norm(full(H_u-H));
display(['the residual W: ', num2str(W_res),' S: ',num2str(S_res),' H ',num2str(H_res)]);
end

[val_H, ind_H]=max(H,[],2);
display('The end of NMTFMG\n');

% f_res=norm(W-(H0)*S*(H1'))+alpha*trace((H0')*Theta0*H0)+beta*trace((H1')*Theta1*H1)+eta1*norm(S(:),1)
% pause;

