function [ R ] = WeightBetweenAttribute(Data)

%   Function:   该函数求出两列间的互信息，并标准化
%   Input:      数据集Data
%   Output:     权重矩阵

[row,col] = size(Data);

for i = 1:col
    for j = i:col
        I = 0;                              %   Mutual Information
        H = 0;                              %   The joint entropy
        Element_i = unique(Data(:,i));      %   取出某一列的元素出现的集合
        Element_j = unique(Data(:,j));      %   取出某一列的元素出现的集合
        
        for Element_i_index = 1:size(Element_i,1)
            for Element_j_index = 1:size(Element_j,1)
                F_i = find(Data(:,i) == Element_i(Element_i_index));
                F_j = find(Data(:,j) == Element_j(Element_j_index));
                P_i = size(F_i,1)/row;     %   计算出在对应列上值等于Element_i(Element_i_index)的元素的个数
                P_j = size(F_j,1)/row;     %   计算出在对应列上值等于Element_i(Element_i_index)的元素的个数
                Temp_i_j = intersect(F_i,F_j);                                       %  在Temp_i的基础上找出值为Element_j(Element_j_index)的元素
                P_i_j = size(Temp_i_j,1)/row;
                if P_i_j == 0           %   没有交集时，互信息为零
                    I = I;
                    H = H;
                else                    %   有交集时，根据互信息公示进行计算
                    I = I + P_i_j * log2(P_i_j/(P_i*P_j));
                    H = H + P_i_j * log2(P_i_j);
                end
            end
        end
        H = -H;
        if H == 0
            R(i,j) = 1;
            R(j,i) = 1;
        else
            R(i,j) = I / H;
            R(j,i) = I / H;
        end
    end
end

[row,col] = size(R);
%   都除以col-1，保证外耦合范围为（0~1）
for i = 1:row
    for j = 1:col
        if i ~=j
            R(i,j) = R(i,j) / (col-1);
        end
    end
end

%   对权重矩阵选取所有列间权重的平均值作为阈值
SumWeight = 0;              %   权重值总和
for i = 1:row
    for j = 1:col
        if i ~=j
            SumWeight = SumWeight + R(i,j);
        end
    end
end
ThresholdValue = SumWeight / (row * row - row);

for i = 1:row
    for j = 1:col
        if R(i,j) < ThresholdValue
            R(i,j) = 0;
            R(j,i) = 0;
        end
    end
end

end