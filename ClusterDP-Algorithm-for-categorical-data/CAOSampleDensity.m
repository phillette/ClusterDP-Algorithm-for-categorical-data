function [ rho ] = CAOSampleDensity( Data )

%   Author: wenjie
%   Data:   2017-8-1
%   功能;   实现分类型数据样本的密度计算

[row,col] = size(Data);

for i = 1:row
    sum = 0;
    for j = 1:col
        persent = size(find(Data(:,j) == Data(i,j)),1) / row;
        sum = sum + persent;
    end
    rho(i) = sum / col;
end

end

