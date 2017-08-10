function [ categoryid ] = WeightAPercentOfMCADO_KMode(Data,K)

% Author Cfy
% Date 2008-02-13
% Modify Date 2008-12-21
% Author WenJie
% Modify Date 2016-09-29
% 参数说明 ：Data 为输入的数据集，InitialCenter是初始聚类中心,K为要求聚类的个数
% Function :实现子文件的K-mode算法，模仿K-means算法,返回样本标签

[n,d] = size(Data);
% 设置categoryid为分类结果的记录，表示向量所属的类
categoryid = zeros(1,n);

% InitialCenters为随机产生的起始中心
InitialCenters = randperm(n);
InitialCenters = InitialCenters(1:K);
% Make this different to get the loop started.
ClusterCenters_i = Data(InitialCenters,[1:d-1]);

while 1 == 1
    %计算每个数据到聚类中心的距离
    for i = 1:n
        dist = Distance_of_Categorical(Data(:,[1:d-1]),i,InitialCenters);% 调用两个对象之间的距离函数
        [m,ind] = min(dist(:,2));                   % 将当前聚类结果存入categoryid中
        categoryid(i) = ind;
    end
    
    ClusterCenters_j = ClusterCenters_i;
    for i = 1:K
        %找到每一类的所有数据，计算他们的Mode，作为下次计算的聚类中心 K-mode算法
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


[row,column] = size(Data);
[a,n] = size(InitialCenters);
dataset = [];
for j = 1:n
    CADO = 0;
    for attribute_index = 1:column
        IaASV = IaASV(Data,Object_i,InitialCenters(j),attribute_index);
        IeASV = IeASV(Data,Object_i,InitialCenters(j),attribute_index);
        CASV = IaASV * IeASV;
        CADO = CADO + CASV;
    end
    dataset = cat(1,dataset, [j,CADO]);
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
    ColumnValue = unique(Data(:,i));      % 得到每一列中不同的属性值
    [Colrow, ColCol] = size(ColumnValue) ;% 得到当前属性值得维数
    ColumnValueSum = 0;
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

function IntraCoupledDissimilarityValue = IaASV(Data,Object_i,Object_j,attribute)

%   Author: wenjie
%   Data:   2017-1-15
%   功能：计算分类型数据集Data中两个对象Object_i,Object_j在属性attribute上的内耦合系数
%   输入参数：数据集Data,对象编号Object_i,对象编号Object_j,属性列attribute
%   输出参数：两个对象Object_i和Object_j在属性列attribute的内耦合系数

global ps;
global pf;

[row,col] = size(Data);

a = size(find(Data(:,attribute) == Data(Object_i,attribute)),1);
b = size(find(Data(:,attribute) == Data(Object_j,attribute)),1);

IntraCoupledSimilarityValue = (a*b)/(a+b+a*b);

if Data(Object_i,attribute) ==  Data(Object_j,attribute)
    weight = ps(attribute) * (a/row) * (b/row);
else
    weight = pf(attribute) * (a/row) * (b/row);
end

IntraCoupledDissimilarityValue = 1 / IntraCoupledSimilarityValue - 1;          %   不相似性
IntraCoupledDissimilarityValue = IntraCoupledDissimilarityValue * weight;

end


function InterCoupledDissimilarityValue = IeASV(Data,Object_i,Object_j,attribute)

%   Author: wenjie
%   Data:   2017-1-15
%   功能：计算分类型数据集Data中两个对象Object_i,Object_j在属性attribute上的相互耦合系数
%   输入参数：数据集Data,对象编号Object_i,对象编号Object_j,属性列attribute
%   输出参数：两个对象Object_i和Object_j在属性列attribute的相互耦合系数
global weight;

[row,col] = size(Data);
%   对象Object_i所在的行与列
[i_row,i_col] = find(Data(:,attribute) == Data(Object_i,attribute));
%   对象Object_j所在的行与列
[j_row,i_col] = find(Data(:,attribute) == Data(Object_j,attribute));

InterCoupledSimilarityValue = 0;
for j = 1:col
    IRSI = 0;
    if j ~= attribute
        %   找到两个对象在非本属性列上的交集
        inter_set = intersect(Data(i_row,j),Data(j_row,j));
        if ~isempty(inter_set)
            for k = 1:size(inter_set,1)
                Object_i_ICP = size(find(Data(i_row,j)==inter_set(k)),1)/size(i_row,1);
                Object_j_ICP = size(find(Data(j_row,j)==inter_set(k)),1)/size(j_row,1);
                IRSI = IRSI + min(Object_i_ICP,Object_j_ICP);
            end
        end
        InterCoupledSimilarityValue = InterCoupledSimilarityValue + weight(j,attribute)*IRSI;     %   相互耦合相似性
    end
end

%   相互耦合不相似性
InterCoupledDissimilarityValue = sum(weight(:,attribute)) - 1 - InterCoupledSimilarityValue;

if abs(InterCoupledDissimilarityValue) < 1 * 10^(-6)
    InterCoupledDissimilarityValue = 0;
end

end
