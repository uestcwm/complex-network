function [I,R]=sire(A,InitialState,beta,mu,time)
%******************************
% A�ڽӾ���
% InitialState��ʼ��Ⱦ״̬
% beta��Ⱦ��
% mu�ָ���
% time����ʱ��
%***************************** 
   Infected=InitialState;recover=[];
   Infected_temp=zeros(size(Infected));
for t=1:1:time %һ������ʱ��time���ݻ�
%%
%��i���׸нڵ�, ���i��һ�����ʽ��д�Ⱦ
    x1=find(Infected==0);
    a1=rand(size(x1));%һ����x1�ߴ���ͬ���������
    b1=beta*(A(x1,:)*Infected);%�Ѿ�����Ⱦ�Ľڵ���һ������ȥ��Ⱦ�����ڵ�
    xx1=setdiff(find(a1<b1),recover);%������A����B�е�Ԫ��  %%�ⲿ�ָ�Ⱦ
    Infected_temp(x1(xx1))=1;%����recover�Ľڵ���һ�����ʱ���Ⱦ
    xx2=setdiff(find(a1>=b1),recover); %%
    Infected_temp(x1(xx2))=0;%�Բ���recover��û����Ⱦ�������ڵ㱣���׸�״̬
%%
%��i��Ⱦ���ڵ�, ���i��һ�����ʽ����Ƴ���recover��
    x2=find(Infected==1); %������Ѿ�����Ⱦ�Ľڵ�
    a2=rand(size(x2));
    xx3=find(a2<mu);%��Ⱦ���ڵ���һ�����ʽ���recover
    recover=[recover;xx3];%���±�recover�ڵ�;
    xx4=find(a2>=mu);
    Infected_temp(x2(xx4))=1;%����δ��recover�Ľڵ��Ա�����Ⱦ����
    Infected=Infected_temp;
    I(t)=sum(Infected); %��¼ÿ��ʱ�䲽��Ⱦ���ڵ�����
    R(t)=sum(recover);
end
end
