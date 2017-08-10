function v = NMI(label, result)

%   Function:   该函数求出两列间的互信息，并标准化
%   Input:      类标签，聚类结果
%   Output:     NMI值

assert(length(label) == length(result));
label = label(:);
result = result(:);

row = length(label);

if length(unique(label)) ~= length(unique(result))
    v = 0;
else
    I = 0;                              %   Mutual Information
    H = 0;                              %   The joint entropy
    HLabel = 0;                         %   Label的熵
    HResult = 0;                        %   Result的熵
    Element_i = unique(label);          %   取出某一列的元素出现的集合
    Element_j = unique(result);         %   取出某一列的元素出现的集合
    
    %   求Label的熵
    for Element_i_index = 1:size(Element_i,1)
        P_i = size(find(label == Element_i(Element_i_index)),1)/row;
        HLabel = HLabel + P_i * log2(P_i + eps);
    end
    HLabel = -HLabel;
    %   求Result的熵
    for Element_j_index = 1:size(Element_j,1)
        P_j = size(find(result == Element_j(Element_j_index)),1)/row;
        HResult = HResult + P_j * log2(P_j + eps);
    end
    HResult = -HResult;
    
    for Element_i_index = 1:size(Element_i,1)
        for Element_j_index = 1:size(Element_j,1)
            F_i = find(label == Element_i(Element_i_index));
            F_j = find(result == Element_j(Element_j_index));
            P_i = size(find(label == Element_i(Element_i_index)),1)/row;        %   计算出在对应列上值等于Element_i(Element_i_index)的元素的个数
            P_j = size(find(result == Element_j(Element_j_index)),1)/row;       %   计算出在对应列上值等于Element_i(Element_i_index)的元素的个数
            Temp_i_j = intersect(F_i,F_j);                                       %  在Temp_i的基础上找出值为Element_j(Element_j_index)的元素
            P_i_j = size(Temp_i_j,1)/row;
            if P_i_j == 0           %   没有交集时，互信息为零
                I = I;
            else                    %   有交集时，根据互信息公示进行计算
                I = I + P_i_j * log2(P_i_j/(P_i*P_j));
            end
        end
    end
    v = 2 * I / (HLabel + HResult);
end
end