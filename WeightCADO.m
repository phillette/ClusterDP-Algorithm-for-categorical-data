function [ DistanceMatrix ] = WeightCADO( Data)

%   Author: wenjie
%   Data:   2017-1-15
%   功能：计算分类型数据集Data中对象间的距离,相比于CADO算法，这里采用Weight替代CADO算法中的1/(col-1)
%   输入参数：分类型数据集Data
%   输出参数：距离矩阵
[row,col] = size(Data);

for i = 1:row
    for j = 1:row
        CADO = 0;
        for attribute_index = 1:col
            IaASV = IaASV(Data,i,j,attribute_index);
            IeASV = IeASV(Data,i,j,attribute_index);
            CASV = IaASV * IeASV;
            CADO = CADO + CASV;
        end
        DistanceMatrix(i,j) = CADO;
        DistanceMatrix(j,i) = CADO;
    end
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
    weight = ps(attribute);
else
    weight = pf(attribute);
end

IntraCoupledDissimilarityValue = 1 / IntraCoupledSimilarityValue - 1;       %   不相似性
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


