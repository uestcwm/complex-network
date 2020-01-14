function [I,R]=sire(A,InitialState,beta,mu,time)
%******************************
% A邻接矩阵
% InitialState初始感染状态
% beta感染率
% mu恢复率
% time传播时间
%***************************** 
   Infected=InitialState;recover=[];
   Infected_temp=zeros(size(Infected));
for t=1:1:time %一共进行时长time的演化
%%
%若i是易感节点, 则对i以一定概率进行传染
    x1=find(Infected==0);
    a1=rand(size(x1));%一个与x1尺寸相同的随机矩阵
    b1=beta*(A(x1,:)*Infected);%已经被感染的节点以一定概率去感染其他节点
    xx1=setdiff(find(a1<b1),recover);%返回在A不在B中的元素  %%这部分感染
    Infected_temp(x1(xx1))=1;%不是recover的节点以一定概率被感染
    xx2=setdiff(find(a1>=b1),recover); %%
    Infected_temp(x1(xx2))=0;%对不是recover且没被感染的其他节点保留易感状态
%%
%若i是染病节点, 则对i以一定概率进行移除（recover）
    x2=find(Infected==1); %如果是已经被感染的节点
    a2=rand(size(x2));
    xx3=find(a2<mu);%对染病节点以一定概率进行recover
    recover=[recover;xx3];%更新被recover节点;
    xx4=find(a2>=mu);
    Infected_temp(x2(xx4))=1;%余下未被recover的节点仍保留感染能力
    Infected=Infected_temp;
    I(t)=sum(Infected); %记录每个时间步的染病节点数量
    R(t)=sum(recover);
end
end
