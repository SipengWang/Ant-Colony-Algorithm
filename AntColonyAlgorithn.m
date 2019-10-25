%��Ⱥ�㷨
clear;
clc;
M = 100;%��������
N = 50;%��������
map(:,1) = randperm(N,N);
map(:,2) = randperm(N,N);
%��������ʼֵ
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
    
    %ΪMֻ�������ѡ��������
    for i = 1:M
        path(i,1) = randperm(N,1);
    end
    
    for i = 1:M %ÿֻ����
        
        allow_reach = [1:N];
        city_index = path(i,1);
        
        for j = 2:N %�������֮���ÿ������
            
            allow_reach(city_index) = 0; %ȥ����һ���Ѿ����ʹ��ĳ���
            
            temp_allow_length = 1; %�õ��ڱ���ѭ����û�б����ʹ��ĳ���
            temp_allow = zeros(N - j + 1, 1);
            for k = 1:N
                if allow_reach(k) ~= 0
                    temp_allow(temp_allow_length) = allow_reach(k);
                    temp_allow_length = temp_allow_length + 1;
                end
            end
            
            %���������Щ���еĸ��ʾ���
            %����һ���Ѿ����ʹ��ĳ��е�����δ���ʳ��еĸ���
            P = zeros(length(temp_allow),1);
            for k = 1:length(temp_allow)
                P(k) = tau(city_index, temp_allow(k))^ALPHA * nt(city_index, temp_allow(k))^BETA;
            end
            
            P = P / sum(P);%���ʹ�һ��.
            P_ = cumsum(P);%�����ۼӣ�������һ����ȡҪ���ʵĳ���
            target = find(P_ >= rand());%�ҵ�P_�����д���rand()��Ԫ�ص��±�
            
            [x, target] = max(P);%
            
            city_index = temp_allow(target(1));%ȡ��һ��
            path(i,j) = city_index;
            
        end
        
    end
    
    len = zeros(M,1);
    %������������ߵ�·���ܳ���
    for i = 1:M
        for j = 1:N - 1
            len(i) = len(i) + D(path(i,j),path(i,j + 1));
        end
        len(i) = len(i) + D(path(i,j + 1),path(i,1));%�������һ�����лص����ľ���
    end
    
    %�������·��
    [min_length, min_index] = min(len);
    best_length(count) = min_length;
    best_path(count,:) = path(min_index,:);
    
    %������Ϣ��
    delta_tau = zeros(N,N);%��Ϣ�ر仯����
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

