function C = NormalizedSC(W, k)

%   Author WenJie
%   Modify Date 2016-09-29
%   该函数实现对聚类的集成，谱聚类算法
%   输入参数Data为类之间距离的矩阵，k表示最后聚类的个数,a代表高斯核函数的参数

[n,m] = size(W);
s = sum(W);
D = full(sparse(1:n, 1:n, s));
E = D^(-1/2)*W*D^(-1/2);
[Q, V] = eigs(E, k);
C = kmeans(Q, k);

end