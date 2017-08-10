function [cl] = Cluster_DP(Data,dist,ClusterNum)

%   doc：       Clustering by fast search and find of density peaks
%   Modify:     2017.2.21
%   Author:     wenjie

global fid;
isShowPicture = 0;          %   是否需要显示聚类图，0表示不显示，1表示显示聚类图


%   将传入的对称的类间矩阵转化为第一列为i,第二列为j，第三列为d（ij）的形式
xx = [];
[row,col] = size(dist);
for i = 1:row
    for j = i+1:col;
        xx = [xx;[i,j,dist(i,j)]];
    end
end

ND = max(xx(:,2));                      %   样本个数
NL = max(xx(:,1));
if NL > ND
    ND = NL;
end
maxd = max(max(dist));              %   密度最大的点，采用max(dij)作为其该点的delta值
nneigh = zeros(1,row);

for bindwidth = 0.1:0.1:0.9
    rho = CateSampleDensity(Data(:,[1:size(Data,2)-1]),bindwidth);
    fprintf(fid,'bindwidth is %.1f \n',bindwidth);
    
    [rho_sorted,ordrho] = sort(rho,'descend');
    delta(ordrho(1)) = -1.;             %   密度最大的点，delta值置为-1
    nneigh(ordrho(1)) = 0;
    
    %   求出每个数据节点的delta值：即与更高密度点的最小距离
    for i = 2:ND
        delta(ordrho(i)) = maxd;
        for j = 1:i-1
            if(dist(ordrho(i),ordrho(j)) < delta(ordrho(i)))
                delta(ordrho(i)) = dist(ordrho(i),ordrho(j));
                %   nneigh记下与该每个数据节点的距离更高密度点的标号
                nneigh(ordrho(i)) = ordrho(j);
            end
        end
    end
    delta(ordrho(1)) = max(delta(:));
    
    if isShowPicture == 1
        scrsz = get(0,'ScreenSize');
        figure('Position',[0 0 scrsz(3) scrsz(4)]);
    end
    
    for i=1:ND
        gamma(i) = rho(i) * delta(i);             %   求出每一个点的rho与delta的乘积,按照从大到小即可得到聚类中心点
    end
    
    if isShowPicture == 1
        %   作图画出Decision Graph图
        subplot(2,1,1);
        plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
        title ('Decision Graph','FontSize',15.0)
        xlabel ('\rho')
        ylabel ('\delta')
    end
    
    %   ordrho记录下rho_sorted在对应于gamma中的索引
    [rho_sorted,ordrho] = sort(gamma,'descend');
    %   ClusterCenterInd记录下gamma中前ClusterNum项作为聚类中心点索引
    ClusterCenterInd = ordrho([1:ClusterNum]);
    
    %   cl标注出该样本点为第几个聚类中心点,-1为非聚类中心点
    for i = 1:ND
        if ismember(i,ClusterCenterInd) == 0        %   当i不是ClusterCenterInd中的元素时
            cl(i) = -1;
        else
            cl(i) = find(i == ClusterCenterInd);    %   当i是ClusterCenterInd中的元素时,记录下是第几个聚类中心
        end
    end
    
    %   记下聚类中心点的样本编号,第i个聚类中心点为源样本中的第几个样本
    for i = 1:ClusterNum
        icl(i) = ClusterCenterInd(i);
    end
    
    %   根据传递规则，找到每一个非聚类中心点样本所在的聚类中心
    while ismember(-1,cl) == 1
        for i = 1:ND
            if cl(ordrho(i)) == -1
                cl(ordrho(i))=cl(nneigh(ordrho(i)));
            end
        end
    end
    
    for i = 1:ClusterNum
        nc = 0;
        nh = 0;
        for j = 1:ND
            %   nc表示第i类聚类结果的样本个数（包括Cluster Core中样本个数和Clister Hole中的样本个数）
            if (cl(j)==i)
                nc = nc + 1;
            end
        end
    end
    
    if isShowPicture == 1
        %   在Decision Graph中画出聚类中心点，颜色区分
        cmap = colormap;
        for i = 1:ClusterNum
            ic = int8((i*64.)/(ClusterNum*1.));
            subplot(2,1,1)
            hold on
            plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        end
    end
    
    if isShowPicture == 1
        subplot(2,1,2)
        Y1 = mdscale(dist, 2, 'criterion','metricsstress');
        plot(Y1(:,1),Y1(:,2),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
        title ('2D Nonclassical multidimensional scaling','FontSize',15.0);
        xlabel ('X');
        ylabel ('Y');
        
        for i = 1:ND
            A(i,1) = 0.;
            A(i,2) = 0.;
        end
        for i = 1:ClusterNum
            nn = 0;
            ic = int8((i*64.)/(ClusterNum*1.));
            for j = 1:ND
            end
            if isShowPicture == 1
                hold on
                plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
            end
        end
    end
    
    [Cluster_DP_AC,Cluster_DP_PR,Cluster_DP_RE,Cluster_DP_CV] = AC_PR_RE(cl,Data(:,size(Data,2)));
    [Cluster_DP_NMI] = NMI(cl,Data(:,size(Data,2)));
    [Cluster_DP_ARI] = AdjustedRandIndex(cl,Data(:,size(Data,2)));
    fprintf(fid,'DP_Average_AC  = %8.4f		DP_Average_NMI = %8.4f      DP_Average_ARI = %8.4f	  	',Cluster_DP_AC,Cluster_DP_NMI,Cluster_DP_ARI);
    fprintf(fid,'DP_Average_PR = %8.4f      DP_Average_RE = %8.4f\n',Cluster_DP_PR,Cluster_DP_RE);
end

end

