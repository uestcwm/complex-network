function [I,R]=sirself(A,InitialState,beta,mu,time)
% A邻接矩阵
% InitialState初始感染状态
% beta感染率
% mu恢复率
% time传播时间
% a=[0 1 0 0 0 
% 1 0 1 0 0 
% 0 1 0 1 1
% 0 0 1 0 0
% 0 0 1 0 0]
I=[];
R=[]; 
Infected=InitialState;recover=[];%开始
for t=1:time
    x=find(Infected==1);
    Li=length(x);%每一次感染集长度
%     if Li==0%如果没有感染节点退出循环
%         break
%     end
    for i=1:Li%这是每次感染集的遍历
       x1=find(A(x(i),:)==1); %遍历到的节点的邻居节点
       Lx=length(x1); %邻居节点的长度
       x2=randi(Lx); %随机取邻居节点的索引
       %% 是否感染:邻居是否感染
       a1=rand;
       if (a1<=beta)&&(Infected(x1(x2))==0)%感染
           Infected(x1(x2))=1;        
       end%感染率小于没被感染 S
       %% 是否恢复：自己是否恢复
        a2=rand;
        if a2<=mu%恢复
            Infected(x(i))=2;
            recover=[recover;Infected(x(i))];
        end               
end 
        I(t)=sum(Infected~=0);%记录每个时间步的染病节点数量
        R(t)=length(recover);
end
end
