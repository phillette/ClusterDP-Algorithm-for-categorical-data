function [ Dist ] = SMS( Data )

%   Author: wenjie
%   Data:   2017-8-6
%   功能;   简单匹配，计算样本距离，返回距离矩阵

[row,col] = size(Data);
Dist = zeros(row,row);

for i = 1:row
    for j = i+1:row
        dist = 0;
        for k = 1:col
            if Data(i,k) ~= Data(j,k)
                dist = dist + 1;
            end
        end
        Dist(i,j) = dist;
        Dist(j,i) = dist;
    end
end


end

