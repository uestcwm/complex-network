function [I,R]=sirself(A,InitialState,beta,mu,time)
% A�ڽӾ���
% InitialState��ʼ��Ⱦ״̬
% beta��Ⱦ��
% mu�ָ���
% time����ʱ��
% a=[0 1 0 0 0 
% 1 0 1 0 0 
% 0 1 0 1 1
% 0 0 1 0 0
% 0 0 1 0 0]
I=[];
R=[]; 
Infected=InitialState;recover=[];%��ʼ
for t=1:time
    x=find(Infected==1);
    Li=length(x);%ÿһ�θ�Ⱦ������
%     if Li==0%���û�и�Ⱦ�ڵ��˳�ѭ��
%         break
%     end
    for i=1:Li%����ÿ�θ�Ⱦ���ı���
       x1=find(A(x(i),:)==1); %�������Ľڵ���ھӽڵ�
       Lx=length(x1); %�ھӽڵ�ĳ���
       x2=randi(Lx); %���ȡ�ھӽڵ������
       %% �Ƿ��Ⱦ:�ھ��Ƿ��Ⱦ
       a1=rand;
       if (a1<=beta)&&(Infected(x1(x2))==0)%��Ⱦ
           Infected(x1(x2))=1;        
       end%��Ⱦ��С��û����Ⱦ S
       %% �Ƿ�ָ����Լ��Ƿ�ָ�
        a2=rand;
        if a2<=mu%�ָ�
            Infected(x(i))=2;
            recover=[recover;Infected(x(i))];
        end               
end 
        I(t)=sum(Infected~=0);%��¼ÿ��ʱ�䲽��Ⱦ���ڵ�����
        R(t)=length(recover);
end
end
