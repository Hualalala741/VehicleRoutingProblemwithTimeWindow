%% ��ջ�������
clear all
clc

%% ��������
load T.mat
load dis.mat
load time.mat
load undi.mat
load xxxx.mat
T=T;
D=dis;
% D=[0,20,30,15;20,0,50,90;30,50,0,40;15,90,40,0]
%% ������м��໥����
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

%% ��ʼ������
m = 1;                              % ��������
alpha = 1;                           % ��Ϣ����Ҫ�̶�����
beta = 5;                            % ����������Ҫ�̶�����
rho = 0.1;                           % ��Ϣ�ػӷ�����
Q = 1;                               % ��ϵ��
Eta = 1./D;                          % ��������
Tau = ones(n,n);                     % ��Ϣ�ؾ���
Table = zeros(m,n);                  % ·����¼��
iter = 1;                            % ����������ֵ
iter_max = 1;                      % ���������� 
Route_best = zeros(iter_max,n);      % �������·��       
Length_best = zeros(iter_max,1);     % �������·���ĳ���  Y
Time_best = zeros(iter_max,1);       % �������·����ʱ��  
Dis_best = zeros(iter_max,1);        % �������·���ĳ���
chi_c_b = zeros(iter_max,1); %�ٵ�����
zao_c_b = zeros(iter_max,1);%�絽����
Length_ave = zeros(iter_max,1);      % ����·����ƽ������  

%% ����Ѱ�����·��
while iter <= iter_max
    % ��������������ϵ�������
      start = zeros(m,1);
      for i = 1:m
          temp = randperm(3)+1;
          temp=[1,temp];
          start(i) = temp(1);%50ֻ����ѡ��
      end
      Table(:,1) = start; %��һ����
      % ������ռ�
      citys_index = 1:n;
      % �������·��ѡ��
      for i = 1:m
          % �������·��ѡ��
         for j = 2:n
             tabu = Table(i,1:(j - 1));           % �ѷ��ʵĳ��м���(���ɱ�)
%              disp('-----')
%              disp(tabu);
             allow_index = ~ismember(citys_index,tabu);
             allow = citys_index(allow_index);  % �����ʵĳ��м���
             P = allow;
             % ������м�ת�Ƹ���
             for k = 1:length(allow)
                 P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % ���̶ķ�ѡ����һ�����ʳ���
             Pc = cumsum(P);     
            target_index = find(Pc >= rand); 
            target = allow(target_index(1));
            Table(i,j) = target;
         end
%         disp(Table);
      end
      % ����������ϵ�·������
      %����ѳ��������ˣ������������
      %�����Ŀ�꺯��
      Length = zeros(m,1);%�����ΪY�����
      Dis=zeros(m,1);
      Time=zeros(m,1);
      Un=zeros(m,1);
      chi_cf=zeros(m,1);
      chi_c=zeros(m,1);
      zao_c=zeros(m,1);
      for i = 1:m
          Route = Table(i,:);%ȡһ��·��
%           disp(Route);
          for j = 1:(n - 1)
              %̫���ˣ�������ֱ�ӰѾ��뵱ʱ��
              %�������Ŀ�꺯��
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
          Length(i) = Length(i) + Un(i)*(Time(i)+chi_cf(i)+Dis(i));%���ϻص���һ��·��
      end
      % �������·�����뼰ƽ������
      if iter == 1
          [min_Length,min_index] = min(Length);%�ҳ�����̵���һ��·��index
          Length_best(iter) = min_Length; %����ÿһ�ε�������õ�һ��· 
          Length_ave(iter) = mean(Length);%ƽ��
          Time_best(iter) = Time(min_index);
          Dis_best(iter) = Dis(min_index);
          chi_c_b(iter)=chi_c(min_index);
          zao_c_b(iter)=zao_c(min_index);
          Route_best(iter,:) = Table(min_index,:);%��õ�·����ֵ
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
      % ������Ϣ��
      Delta_Tau = zeros(n,n);
      % ������ϼ���
      for i = 1:m
          % ������м���
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-rho) * Tau + Delta_Tau;
    % ����������1�����·����¼��
    iter = iter + 1;
    Table = zeros(m,n);
end

%% �����ʾ
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp('ʱ��')
chi_c_b(index)
zao_c_b(index)
Time_best(index)
disp('����')
Dis_best(index)
disp(['��̾���:' num2str(Shortest_Length)]);
disp(['���·��:' num2str([Shortest_Route Shortest_Route(1)])]);

%% ��ͼ
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       ���');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       �յ�');
xlabel('����')
ylabel('γ��')
title(['��̾���·��'])
figure(2)
plot(1:iter_max,Length_best,1:iter_max,Length_ave,'r:')
legend('����Ŀ�꺯��ֵ','ƽ��ֵ')
xlabel('��������')
ylabel('ֵ')
title('��������ֵ��ƽ��ֵ�Ա�')
