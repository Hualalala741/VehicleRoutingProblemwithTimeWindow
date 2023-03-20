%% 清空环境变量
clear all
clc

%% 导入数据
load T.mat
load dis.mat
load time.mat
load undi.mat
load xxxx.mat
T=T;
D=dis;
% D=[0,20,30,15;20,0,50,90;30,50,0,40;15,90,40,0]
%% 计算城市间相互距离
n=50
% n = size(citys,1);
% D = zeros(n,n);
% for i = 1:n
%     for j = 1:n
%         if i ~= j
%             D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
%         else
%             D(i,j) = 1e-4;      
%         end
%     end    
% end
% D=dist

%% 初始化参数
m = 1;                              % 蚂蚁数量
alpha = 1;                           % 信息素重要程度因子
beta = 5;                            % 启发函数重要程度因子
rho = 0.1;                           % 信息素挥发因子
Q = 1;                               % 常系数
Eta = 1./D;                          % 启发函数
Tau = ones(n,n);                     % 信息素矩阵
Table = zeros(m,n);                  % 路径记录表
iter = 1;                            % 迭代次数初值
iter_max = 1;                      % 最大迭代次数 
Route_best = zeros(iter_max,n);      % 各代最佳路径       
Length_best = zeros(iter_max,1);     % 各代最佳路径的长度  Y
Time_best = zeros(iter_max,1);       % 各代最佳路径的时间  
Dis_best = zeros(iter_max,1);        % 各代最佳路径的长度
chi_c_b = zeros(iter_max,1); %迟到次数
zao_c_b = zeros(iter_max,1);%早到次数
Length_ave = zeros(iter_max,1);      % 各代路径的平均长度  

%% 迭代寻找最佳路径
while iter <= iter_max
    % 随机产生各个蚂蚁的起点城市
      start = zeros(m,1);
      for i = 1:m
          temp = randperm(3)+1;
          temp=[1,temp];
          start(i) = temp(1);%50只蚂蚁选择
      end
      Table(:,1) = start; %第一个点
      % 构建解空间
      citys_index = 1:n;
      % 逐个蚂蚁路径选择
      for i = 1:m
          % 逐个城市路径选择
         for j = 2:n
             tabu = Table(i,1:(j - 1));           % 已访问的城市集合(禁忌表)
%              disp('-----')
%              disp(tabu);
             allow_index = ~ismember(citys_index,tabu);
             allow = citys_index(allow_index);  % 待访问的城市集合
             P = allow;
             % 计算城市间转移概率
             for k = 1:length(allow)
                 P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % 轮盘赌法选择下一个访问城市
             Pc = cumsum(P);     
            target_index = find(Pc >= rand); 
            target = allow(target_index(1));
            Table(i,j) = target;
         end
%         disp(Table);
      end
      % 计算各个蚂蚁的路径距离
      %上面把城市跑完了，现在来算最短
      %这个是目标函数
      Length = zeros(m,1);%这个作为Y矩阵吧
      Dis=zeros(m,1);
      Time=zeros(m,1);
      Un=zeros(m,1);
      chi_cf=zeros(m,1);
      chi_c=zeros(m,1);
      zao_c=zeros(m,1);
      for i = 1:m
          Route = Table(i,:);%取一条路径
%           disp(Route);
          for j = 1:(n - 1)
              %太难了，我现在直接把距离当时间
              %这个就是目标函数
%               Length(i) = Length(i) + D(Route(j),Route(j + 1));%Y
              Dis(i)=Dis(i)+D(Route(j),Route(j+1));%p1
              Time(i)=Time(i)+time(Route(j),Route(j+1));%p2
              Un(i)=Un(i)+undi(Route(j),Route(j+1));%p4
              %p3
              if Time(i)<T(Route(j+1),1)
                  Time(i)=Time(i)+T(Route(j+1),1)-Time(i);
                  zao_c(i)=zao_c(i)+1;
              elseif Time(i)>T(Route(j+1),2)
                  chi_cf(i)=chi_cf(i)+2*(Time(i)-T(Route(j+1),2));
                  chi_c(i)=chi_c(i)+1;
              end

          end
%           Length(i) = Length(i) + D(Route(n),Route(1));
          Dis(i)=Dis(i)+D(Route(n),Route(1));
          Time(i)=Time(i)+time(Route(n),Route(1));
          Length(i) = Length(i) + Un(i)*(Time(i)+chi_cf(i)+Dis(i));%加上回到第一条路的
      end
      % 计算最短路径距离及平均距离
      if iter == 1
          [min_Length,min_index] = min(Length);%找长度最短的那一条路的index
          Length_best(iter) = min_Length; %储存每一次迭代里最好的一条路 
          Length_ave(iter) = mean(Length);%平均
          Time_best(iter) = Time(min_index);
          Dis_best(iter) = Dis(min_index);
          chi_c_b(iter)=chi_c(min_index);
          zao_c_b(iter)=zao_c(min_index);
          Route_best(iter,:) = Table(min_index,:);%最好的路和他值
      else
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);
          Length_ave(iter) = mean(Length);
          Time_best(iter) = Time(min_index);
          Dis_best(iter) = Dis(min_index);
          chi_c_b(iter)=chi_c(min_index);
          zao_c_b(iter)=zao_c(min_index);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
      % 更新信息素
      Delta_Tau = zeros(n,n);
      % 逐个蚂蚁计算
      for i = 1:m
          % 逐个城市计算
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-rho) * Tau + Delta_Tau;
    % 迭代次数加1，清空路径记录表
    iter = iter + 1;
    Table = zeros(m,n);
end

%% 结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp('时间')
chi_c_b(index)
zao_c_b(index)
Time_best(index)
disp('距离')
Dis_best(index)
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);

%% 绘图
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       起点');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       终点');
xlabel('经度')
ylabel('纬度')
title(['最短距离路径'])
figure(2)
plot(1:iter_max,Length_best,1:iter_max,Length_ave,'r:')
legend('最优目标函数值','平均值')
xlabel('迭代次数')
ylabel('值')
title('各代最优值与平均值对比')
