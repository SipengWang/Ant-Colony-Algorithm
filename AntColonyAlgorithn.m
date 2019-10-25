%蚁群算法
clear;
clc;
M = 100;%蚂蚁数量
N = 50;%城市数量
map(:,1) = randperm(N,N);
map(:,2) = randperm(N,N);
%参量赋初始值
L = 200;
ALPHA = 5;
BETA = 12;
RHO = 0.1;
tau = ones(N,N);
len = zeros(M,1);


for i = 1:N
    for j = 1:i
        
        if i == j
            D(i,j) = 10^(-5);
        else
            D(i,j) = norm(map(i,:) - map(j,:));
        end
        
    end
end

D = D + D';

nt = 1./D;


for count = 1:L
    
    %为M只蚂蚁随机选择起点城市
    for i = 1:M
        path(i,1) = randperm(N,1);
    end
    
    for i = 1:M %每只蚂蚁
        
        allow_reach = [1:N];
        city_index = path(i,1);
        
        for j = 2:N %除了起点之外的每个城市
            
            allow_reach(city_index) = 0; %去除上一个已经访问过的城市
            
            temp_allow_length = 1; %得到在本次循环中没有被访问过的城市
            temp_allow = zeros(N - j + 1, 1);
            for k = 1:N
                if allow_reach(k) ~= 0
                    temp_allow(temp_allow_length) = allow_reach(k);
                    temp_allow_length = temp_allow_length + 1;
                end
            end
            
            %计算访问这些城市的概率矩阵
            %从上一个已经访问过的城市到其他未访问城市的概率
            P = zeros(length(temp_allow),1);
            for k = 1:length(temp_allow)
                P(k) = tau(city_index, temp_allow(k))^ALPHA * nt(city_index, temp_allow(k))^BETA;
            end
            
            P = P / sum(P);%概率归一化.
            P_ = cumsum(P);%概率累加，方便下一步抽取要访问的城市
            target = find(P_ >= rand());%找到P_中所有大于rand()的元素的下标
            
            [x, target] = max(P);%
            
            city_index = temp_allow(target(1));%取第一个
            path(i,j) = city_index;
            
        end
        
    end
    
    len = zeros(M,1);
    %计算各个蚂蚁走的路径总长度
    for i = 1:M
        for j = 1:N - 1
            len(i) = len(i) + D(path(i,j),path(i,j + 1));
        end
        len(i) = len(i) + D(path(i,j + 1),path(i,1));%加上最后一个城市回到起点的距离
    end
    
    %计算最佳路径
    [min_length, min_index] = min(len);
    best_length(count) = min_length;
    best_path(count,:) = path(min_index,:);
    
    %更新信息素
    delta_tau = zeros(N,N);%信息素变化矩阵
    for i = 1:M
        for j = 1:N - 1
            delta_tau(path(i,j),path(i,j + 1)) =  delta_tau(path(i,j),path(i,j + 1)) + 1/len(i);
        end
        delta_tau(path(i,1),path(i,N)) = delta_tau(path(i,1),path(i,N)) + 1/len(i);
    end
    tau = (1 - RHO) * tau + delta_tau;

    
end

[length, index] = min(best_length);

plot(map(best_path(index,:)',1),map(best_path(index,:)',2), '*-')

