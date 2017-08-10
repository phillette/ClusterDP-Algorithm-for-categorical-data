function RI = RandIndex(X,Y)

%
%	假设一个集合中有N篇文章,一个集合中有N(N-1)/2个集合对
%   TP：同一类的文章被分到同一个簇
%   TN：不同类的文章被分到不同簇
%   FP：不同类的文章被分到同一个簇
%   FN：同一类的文章被分到不同簇
% 
%   Rand Index度量的正确的百分比 RI = （TP+TN）/（TP+FP+FN+TN）

r = length(X);
xk = unique(X);
xk_len = length(xk);
yk = unique(Y);
yk_len = length(yk);
Nx = 0;
for i = 1:xk_len
    for j = 1:yk_len
        A = find(X==xk(i));
        B = find(Y==yk(j));
        C = intersect(A,B);
        N(i,j) = length(C);
        Nx = Nx + N(i,j)*(N(i,j)-1)/2;
    end
end
a = zeros(1,xk_len);
Ax = 0;
for i = 1:xk_len
    a(i) = sum(N(i,:));
    Ax = Ax+a(i)*(a(i)-1)/2;
end
Bx = 0;
b = zeros(1,yk_len);
for i = 1:yk_len
    b(i) = sum(N(:,i));
    Bx = Bx+b(i)*(b(i)-1)/2;
end
RI = (Nx-Ax*Bx/(r*(r-1)/2))/((Ax+Bx)/2-Ax*Bx/(r*(r-1)/2));
end