function val=tri_dense_eval(A1, A2, W12,cluster1, cluster2)

average=(sum(sum(A1))+sum(sum(A2))+sum(sum(W12)))/(size(A1,1)+size(A2,1));

dense=(sum(sum(A1(cluster1,cluster1)))+sum(sum(A2(cluster2,cluster2)))+sum(sum(W12(cluster1, cluster2))))/(length(cluster1)+length(cluster2));

val=dense/average;