function C = SpectralClustering(Data, k, a)  

%   Author WenJie
%   Modify Date 2016-09-29
%   该函数实现对聚类的集成，谱聚类算法 
%   输入参数Data为类之间距离的矩阵，k表示最后聚类的个数,a代表高斯核函数的参数

d = pdist(Data); 
d2 = squareform(d);
d3 = d2.^2;
W(:,:) = exp(-d3(:,:)/(2*a^2));   % 用高斯核函数（也称径向基函数核）计算相似度，距离越大，代表其相似度越小
[n,m] = size(W);  
s = sum(W);
D = full(sparse(1:n,1:n,s));
E = D^(-1/2)*W*D^(-1/2);           % 归一化拉普兰斯矩阵 
[X,B] = eig(E);                    % 求出归一化之后的特征值和特征向量
[Q,V] = eigs(E,k);                 % 选的是前K个最大特征值对应的特征向量
C = kmeans(Q,k);                   % 对特征矩阵做K-Means做聚类


end 

