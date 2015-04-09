clear;

A_i=[1;1;1;2;3;4;5;6;5;7;4;7;9;5;9;11;12;12;13;14;13;14;13];
A_j=[2;3;4;4;4;5;6;7;8;8;9;10;11;12;12;13;13;14;14;15;16;16;17];
A_v=[4;4;3;3;4;1;1;1;1;1;1;1;4;1.5;1;1;4;4;3;1;3.5;2;2];

A_i_d=(1:17)';
A_j_d=(1:17)';
A_v_d=3.5*ones(17,1);




A1_i=[1;1;1;2;3;4;5;6;5;7;4;7;9;5;9;11;12;12;13;14;13;14;13;7;7;18;18;19];
A1_j=[2;3;4;4;4;5;6;7;8;8;9;10;11;12;12;13;13;14;14;15;16;16;17;18;20;20;19;20];
A1_v=[4;4;3;3;4;1;1;1;1;1;1;1;4;1.5;1;1;4;4;3;1;3.5;2;2;4;4;4;3;3.5];

A1_i_d=(1:20)';
A1_j_d=(1:20)';
A1_v_d=3.5*ones(20,1);

% A_i=[A_i;A_i_d];
% A_j=[A_j;A_j_d];
% A_v=[A_v;A_v_d];
% A1_i=[A1_i;A1_i_d];
% A1_j=[A1_j;A1_j_d];
% A1_v=[A1_v;A1_v_d];

A0=sparse(A_i,A_j,A_v,17,17);
A1=sparse(A1_i,A1_j,A1_v,20,20);
A0=A0+A0';
A1=A1+A1';

%%1
% W_i=[1;2;3;4;5;1;2];
% W_j=[12;13;14;16;11;3;4];
% W_v=[4;4;4;3;4;1;1];

%%2
% W_i=[1;2;3;4;5;1;2;12;14;13;16];
% W_j=[12;13;14;16;11;3;4;7;18;20;19];
% W_v=[4;4;4;3;4;1;1;1;1;1;1];

%%3
W_i=[1;2;3;4;5;1;2;12;14;13;16];
W_j=[12;13;14;16;11;3;4;7;18;20;19];
W_v=[4;4;4;3;4;1;1;4;4;4;4];

%%synthetic graph 4
% W_i=[1;2;3;4;5;1;2;12;14;13;16;12;13;14;16];
% W_j=[12;13;14;16;11;3;4;7;18;20;19;1;3;4;2];
% W_v=[4;4;4;3;4;1;1;4;4;4;4;1;1;1;1];

% W_i=[1;2;3;4;5;1;2;12;14;13;16;12;13;14;16];
% W_j=[12;13;14;16;11;3;4;7;18;20;19;1;3;4;2];
% W_v=[4;4;4;3;4;1;1;4;4;4;4;4;4;4;4];




W=sparse(W_i,W_j,W_v,17,20);


k0=4; 
k1=4;  

method='RCC';

switch method
    case 'DRCC',
        DRCC(W, A0, A1, k0, k1);
    case 'DNMTF'
        DNMTF(W, A0, A1, k0, k1);
    case 'NMTFOC' 
        NMTFOC(W, A0, A1, k0, k1);
    case 'RCC'
        [F, G H]=RCC(W, A0, A1, k0, k1);
        [val_H0, ind_H0]=max(F,[],2);
        [val_H1, ind_H1]=max(G,[],2);
        
    otherwise,
        ;
end
        
