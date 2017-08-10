function [ps,pf] = WeightInAttribute(Data)

%   Function:   该函数求出每一列熵作为内耦合权重值
%   Input:      数据集Data
%   Output:     权重矩阵

ps = ProbabilitySameValue( Data );
pf = ProbabilityDifferent( Data );

end

function [ ps ] = ProbabilitySameValue( Data )

%   Function:   该函数求出某一列中取相同值的概率
%   Input:      数据集Data
%   Output:     每一列取相同值的概率标量

[row,col] = size(Data);
for i = 1:col
    probability = 0;
    mr = unique(Data(:,i));     %   取出某一列的元素出现的集合
    for j = 1:size(mr,1)
        EqualXiNum = size(find(Data(:,i) == mr(j)),1);     %   计算出在对应列上值等于Xi的元素的个数
        probability = probability + (EqualXiNum / row) * ((EqualXiNum-1) / (row-1));
    end
    ps(i) = probability;
end

end

function [ pf ] = ProbabilityDifferent( Data )

%   Function:   该函数求出某一列中取不同值的概率
%   Input:      数据集Data
%   Output:     每一列取不同值的概率标量

pf = 1.- ProbabilitySameValue(Data);

end