% function spreadingability=spread(A,beta,mu)
clc;clear all;
A=load ('arpa.mat');
A=struct2cell(A);A=cell2mat(A);beta=0.4;mu=0.5;

n=100;%循环量，蒙特卡洛次数
time=20;%传播时间
time_xianshi=40;
for g=1:25%前top-k节点循环
InitialState1=zeros(length(A),1);
% InitialState1(14)=1;InitialState1(3)=1;        
% InitialState1(6)=1;InitialState1(12)=1;
InitialState1(index(1:g))=1;%可以认为指定前几个节点
    for N=1:n%随机次数
        I(N,:)=sire(A,InitialState1,beta,mu,time);
    end
%     spreadingability(i,1)=mean(mean(I));%节点i的传播能力
% end
I_mean=mean(I)/length(A);
I_mean_xianshi(g)=I_mean(time_xianshi);
%% %2DC度中心性方法
    InitialState2=zeros(length(A),1);
    InitialState2(z1(1:g))=1;
%     InitialState2(14)=1;InitialState2(3)=1;
%     InitialState2(6)=1;InitialState2(2)=1;
for N=1:n
    I2(N,:)=sire(A,InitialState2,beta,mu,time);
end
I2_mean=mean(I2)/length(A);
I2_mean_xianshi(g)=I2_mean(time_xianshi);
%% %3BC介数中心性方法
    InitialState3=zeros(length(A),1);
    InitialState3(z3(1:g))=1;
%     InitialState3(12)=1;InitialState3(3)=1;
%     InitialState3(6)=1;InitialState3(19)=1;
for N=1:n
    I3(N,:)=sire(A,InitialState3,beta,mu,time);
end
I3_mean=mean(I3)/length(A);
I3_mean_xianshi(g)=I3_mean(time_xianshi);
%% 4CC接近度中心性方法
InitialState4=zeros(length(A),1);
InitialState4(z2(1:g))=1;
% InitialState4(19)=1;InitialState4(3)=1;
% InitialState4(18)=1;InitialState4(12)=1;
for N=1:n%随机次数
    I4(N,:)=sire(A,InitialState4,beta,mu,time);
end
I4_mean=mean(I4)/length(A);
I4_mean_xianshi(g)=I4_mean(time_xianshi);
%% 5SH结构洞方法
    InitialState5=zeros(length(A),1);
    InitialState5(zSH(1:g))=1;
%     InitialState5(3)=1;InitialState5(14)=1;
%     InitialState5(6)=1;InitialState5(12)=1;
for N=1:n%随机次数
    I5(N,:)=sire(A,InitialState5,beta,mu,time);
end
I5_mean=mean(I5)/length(A);
I5_mean_xianshi(g)=I5_mean(time_xianshi);
end
%% 第一种仿真：感染率随时间
% t=1:time;
% plot(t,I_mean,'o-',t,I2_mean,'g-d',t,I3_mean,'r-s',t,I4_mean,'m+-',t,I5_mean,'k*-');
% xlabel('感染时间');ylabel('感染百分比%');
% hleg2=legend('ISH','DC','BC','CC','SH');
% set(hleg2, 'Location', 'SouthEast');
%% 第2种：感染率随初始感染节点
g=1:25;
plot(g,I_mean_xianshi,'o-',g,I2_mean_xianshi,'g-d',...
  g,I3_mean_xianshi,'r-s',g,I4_mean_xianshi,'m+-',g,I5_mean_xianshi,'k*-');
xlabel('初始感染节点');ylabel('感染百分比%');
hleg2=legend('ISH','DC','BC','CC','SH');
set(hleg2, 'Location', 'SouthEast');