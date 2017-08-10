function [ categoryid ] = NDM_ForCD_K_Mode(Data,K)

% Author WenJie
% Modify Date 2017-04-11
% 参数说明 ：Data 为输入的数据集，K为要求聚类的个数
% Function :实现子文件的K-mode算法，模仿K-means算法,返回样本标签

[n,d] = size(Data);
% 设置categoryid为分类结果的记录，表示向量所属的类
categoryid = zeros(1,n);

% InitialCenters为随机产生的起始中心
InitialCenters = randperm(n);
InitialCenters = InitialCenters(1:K);
%   Make this different to get the loop started.
ClusterCenters_i = Data(InitialCenters,[1:d-1]);

while 1 == 1
    %   计算每个数据到聚类中心的距离
    for i = 1:n
        Dist = Distance_of_Categorical(Data(:,[1:d-1]),i,InitialCenters);   %  调用两个对象之间的距离函数
        [m,ind] = min(Dist(:,2));                                           %   将当前聚类结果存入categoryid中
        categoryid(i) = ind;
    end
    
    ClusterCenters_j = ClusterCenters_i;
    for i = 1:K
        %   找到每一类的所有数据，计算他们的Mode，作为下次计算的聚类中心K-mode算法
        ind = find(categoryid == i);
        if ~isempty(ind)
            ClusterCenters_i(i,:) = Find_Mode(Data(ind,[1:d-1]));
        end
    end
    
    if ClusterCenters_i == ClusterCenters_j
        %     % 统计每一类的数据个数
        %     for i = 1:K
        %         ind = find(categoryid == i+VectorIndex);
        %         CategoryLength = [CategoryLength,length(ind)];
        %     end
        break
    end
end

end

function [dataset]=Distance_of_Categorical(Data,Object_i,InitialCenters)

% Author Fuyuan Cao
% Data : 2008-02-13
% Function :
% 计算两个具有相同的行列的矩阵中，相同行的距离

[row,col] = size(Data);
[a,n] = size(InitialCenters);
dataset = [];
global weight;
global ps;
global pf;

for k = 1:n
    Dist = 0;
    d = 0;
    w = 0;
    for i = 1:col
        if Data(Object_i,i) ~= Data(InitialCenters(k),i)
            for j = 1 : col
                F_i_i = find(Data(:,i) == Data(Object_i,i));
                F_i_j = find(Data(:,j) == Data(Object_i,j));
                Temp_i = intersect(F_i_i,F_i_j);
                P_i = size(Temp_i,1)/row;
                P_i_No = (size(Temp_i,1) - 1)/(row - 1);
                F_j_i = find(Data(:,i) == Data(InitialCenters(k),i));
                F_j_j = find(Data(:,j) == Data(InitialCenters(k),j));
                Temp_j = intersect(F_j_i,F_j_j);
                P_j = size(Temp_j,1)/row;
                P_j_No = (size(Temp_j,1) - 1)/(row - 1);
                d = d + weight(i,j) * (P_i * P_i_No + P_j * P_j_No);
            end
            wr = ps(i);
        else
            for j = 1 : col
                F_i_i = find(Data(:,i) == Data(Object_i,i));
                F_i_j = find(Data(:,j) == Data(Object_i,j));
                Temp_i = intersect(F_i_i,F_i_j);
                P_i = size(Temp_i,1)/row;
                P_i_No = (size(Temp_i,1) - 1)/(row - 1);
                F_j_i = find(Data(:,i) == Data(InitialCenters(k),i));
                F_j_j = find(Data(:,j) == Data(InitialCenters(k),j));
                Temp_j = intersect(F_j_i,F_j_j);
                P_j = size(Temp_j,1)/row;
                P_j_No = (size(Temp_j,1) - 1)/(row - 1);
                d = d + weight(i,j) * (P_i * P_i_No + P_j * P_j_No);
            end
            wr = pf(i);
        end
        w = w + wr;
        Dist = Dist + wr*d;
    end
    Dist = Dist / w;
    dataset = cat(1,dataset, [k,Dist]);
end
end


function [Mode] = Find_Mode(Data)

% Author Cfy
% Data 2008-02-13
% 参数说明 ：Data 为输入的数据集
% Function : 得到一个数据集合的mode

[row,column] = size(Data);
Mode = [];
for i=1:column
    ColumnValue = unique(Data(:,i));        % 得到每一列中不同的属性值
    [Colrow, ColCol] = size(ColumnValue) ;  % 得到当前属性值得维数
    ColumnValueSum = 0;                     % 表示第j个属性值在第i列出现的次数
    dataset = [];
    for j = 1:Colrow
        [Rowresult,Colresult] = size(find(Data(:,i)==ColumnValue(j,1)));
        %   ColumnValueSum = Rowresult;           % 表示第j个属性值在第i列出现的次数
        %   ColumnValue(j,2) = ColumnValueSum     % 得到每一个属性值出现的次数
        dataset = cat(1,dataset,[ColumnValue(j,1),Rowresult]);
    end
    [Maxrow,MaxValue] = max(dataset(:,2));
    Mode(1,i) = ColumnValue(MaxValue,1);
end
end

function [ delta ] = delta( Xi , Xj )

%   Function:   该函数求出同一列属性中两个值间的权重值
%   Input:      样本编号Xi,Xj
%   Output:     两个样本之间的delta值

if Xi ~= Xj
    delta = 1;
else
    delta = 0;
end

end
