clc;
clear all;
close all;
tic
[a,Du,txt_line_number] = txt(  );%email.txt转换为mat邻接矩阵
N=size(a,1);
%% 1
DC=Du/(N-1);
%% 
BB=sparse(a);
d=graphallshortestpaths(BB,'Directed',0); %求所有顶点对之间的最短距离
%% 2k-shell
ks=mycoreness(a);%ascend
%% 3Cnc
Cnc=zeros(1,N);
for i=1:N
    for j=1:N
        if a(i,j)~=0
        Cnc(i)=Cnc(i)+ks(j);
        end
    end
end 
Cncopt=zeros(1,N);
for i=1:N
    for j=1:N
        if a(i,j)~=0
        Cncopt(i)=Cncopt(i)+Cnc(j);
        end
    end
end 
%% 5
w=zeros(N,N);
for i=1:N
    for j=1:N
        w(i,j)=Du(i)*Du(j);
    end
end
for i=1:N
    for j=1:N
        if a(i,j)==0
            w(i,j)=0;
        end
    end
end
w_mean=sum(sum(w))/(2*txt_line_number);%w_mean=(sum(sum(w))-sum(sum(tril(w))))/txt_line_number
C_fai=zeros(1,N);
for i=1:N
    for j=1:N
        if a(i,j)~=0
            C_fai(i)=C_fai(i)+w(i,j)*Du(j)./w_mean;
        end
    end
    C_fai(i)=C_fai(i)+Du(i);
end    
%% 6 ksk
ksk=ks_k(a,N,d,ks);
%% 7 MDD
[result_km,M_mdd,km,res_index,append]=mdd(a,N);
%% 8 ksh
ksh = kshreview_4(a,N,Du,ks);
%% 9IKS
[result_iks,index,result9,ze,e] = IKS(a,N,Du,ks);
%% 10
mcde=MCDE( a,Du,ks,N);
%% 
[result1,z1] = sort(DC,'descend');
[result2,z2] = sort(ks,'descend');
[result3Cnc,zCnc] = sort(Cnc,'descend');%
[result4Cncopt,zCncopt] = sort(Cncopt,'descend');%
[result5C_fai,zC_fai] = sort(C_fai,'descend'); %
[result6ksk,zksk] = sort(ksk,'ascend');%
[result7km,zkm] = sort(result_km,'descend');%
[result8ksh,zksh] = sort(ksh,'descend');%
[result10mcde,zmcde] = sort(mcde,'descend');%
%% 开始验证性能
beta=0.02;mu=0.01; 
n=700;%循环量，蒙特卡洛次数
time=600;
time_xianshi=600;
nlen=60;% nlen=floor(N*0.2);%
data1=cell(1,n);data2=cell(1,n);data3=cell(1,n);data4=cell(1,n);data5=cell(1,n);
data6=cell(1,n);data7=cell(1,n);data8=cell(1,n);data9=cell(1,n);
I1data=0;I2data=0;I3data=0;I4data=0;I5data=0;I6data=0;I7data=0;I8data=0;I9data=0;
for g=1:floor(nlen/20):nlen%top-k
%% 1DC
    InitialState1=zeros(length(a),1);
    InitialState1(z1(1:g))=1;
    for i=1:n
       [I,R]=sirself(a,InitialState1,beta,mu,time);
       data1{i}=I;%data11{i}=R;
       I1data=I1data+data1{1,i}(1,time_xianshi);
    end
    I1data=I1data/n;
    I_mean_xianshi(floor(g/floor(nlen/20))+1)=I1data/N;%R_mean_xianshi(g)=dat11/(n*length(a));
%     I_mean_xianshi(floor(g/2)+1)=I1data/N;
disp(g);
%% 2KSHELL
    InitialState2=zeros(length(a),1);
    InitialState2(z2(1:g))=1;
    for i=1:n
        [I2,R2]=sirself(a,InitialState2,beta,mu,time);
        data2{i}=I2;% data44{i}=R4;       
        I2data=I2data+data2{1,i}(1,time_xianshi);
    end
    I2data=I2data/n;%R6data=R6data/n;    
    I2_mean_xianshi(floor(g/floor(nlen/20))+1)=I2data/N;%R4_mean_xianshi(g)=R6data/N;
%     I2_mean_xianshi(floor(g/2)+1)=I2data/N;
%% 3Cnc
    InitialState3=zeros(length(a),1);
    InitialState3(zCnc(1:g))=1;
    for i=1:n
        [I3,R3]=sirself(a,InitialState3,beta,mu,time);
        data3{i}=I3; %data66{i}=R6;    
        I3data=I3data+data3{1,i}(1,time_xianshi);
    end
    I3data=I3data/n;%R6data=R6data/n;    
    I3_mean_xianshi(floor(g/floor(nlen/20))+1)=I3data/N;%R6_mean_xianshi(g)=R6data/N;
%     I3_mean_xianshi(floor(g/2)+1)=I3data/N;
%% 4Cnc+
    InitialState4=zeros(length(a),1);
    InitialState4(zCncopt(1:g))=1;
    for i=1:n
        [I4,R4]=sirself(a,InitialState4,beta,mu,time);
        data4{i}=I4; %data77{i}=R7;     
        I4data=I4data+data4{1,i}(1,time_xianshi);
    end
    I4data=I4data/n;%R7data=R7data/n;    
    I4_mean_xianshi(floor(g/floor(nlen/20))+1)=I4data/N;%R7_mean_xianshi(g)=R7data/N;
%     I4_mean_xianshi(floor(g/2)+1)=I4data/N;
%% 5C_fai
    InitialState5=zeros(length(a),1);
    InitialState5(zC_fai(1:g))=1;
    for i=1:n
        [I5,R5]=sirself(a,InitialState5,beta,mu,time);
        data5{i}=I5;      
        I5data=I5data+data5{1,i}(1,time_xianshi);
    end
    I5data=I5data/n;    
    I5_mean_xianshi(floor(g/floor(nlen/20))+1)=I5data/N;
%     I5_mean_xianshi(floor(g/2)+1)=I5data/N;
%% 6 ksk
    InitialState6=zeros(length(a),1);
    InitialState6(zksk(1:g))=1;
    for i=1:n
        [I6,R6]=sirself(a,InitialState6,beta,mu,time);
        data6{i}=I6;% data88{i}=R8;   
        I6data=I6data+data6{1,i}(1,time_xianshi);
    end
    I6data=I6data/n;%R8data=R8data/n;    
    I6_mean_xianshi(floor(g/floor(nlen/20))+1)=I6data/N;%R8_mean_xianshi(g)=R8data/N;
%     I6_mean_xianshi(floor(g/2)+1)=I6data/N;
%% 7 km
    InitialState7=zeros(length(a),1);
    InitialState7(zkm(1:g))=1;
    for i=1:n
        [I7,R7]=sirself(a,InitialState7,beta,mu,time);
        data7{i}=I7; %data99{i}=R9;      
        I7data=I7data+data7{1,i}(1,time_xianshi);
    end
    I7data=I7data/n;%R9data=R9data/n;    
    I7_mean_xianshi(floor(g/floor(nlen/20))+1)=I7data/N;%R9_mean_xianshi(g)=R9data/N; 
%     I7_mean_xianshi(floor(g/2)+1)=I7data/N;
%% 8 ksh
    InitialState8=zeros(length(a),1);
    InitialState8(zmcde(1:g))=1;
    for i=1:n
        [I8,R8]=sirself(a,InitialState8,beta,mu,time);
        data8{i}=I8; %data1010{i}=R10;  
        I8data=I8data+data8{1,i}(1,time_xianshi);
    end
    I8data=I8data/n;%R10data=R10data/n;    
    I8_mean_xianshi(floor(g/floor(nlen/20))+1)=I8data/N;%R10_mean_xianshi(g)=R10data/N;
%     I8_mean_xianshi(floor(g/2)+1)=I8data/N;
%% 9IKS
    InitialState9=zeros(length(a),1);
    InitialState9(index(1:g))=1;
    for i=1:n
        [I9,R9]=sirself(a,InitialState9,beta,mu,time);
        data9{i}=I9; %data99{i}=R9;
        I9data=I9data+data9{1,i}(1,time_xianshi);
    end
    I9data=I9data/n;%R9data=R9data/n;    
    I9_mean_xianshi(floor(g/floor(nlen/20))+1)=I9data/N;%R9_mean_xianshi(g)=R9data/N;
%     I9_mean_xianshi(floor(g/2)+1)=I9data/N;
 end
%% 
for i=1:n
    m_1(i)=length(data1{i});
end
m=max(m_1);
%% 
for i=1:n
    matrix1{i}=[data1{i}];matrix2{i}=[data2{i}];matrix3{i}=[data3{i}];matrix4{i}=[data4{i}];matrix5{i}=[data5{i}];
    matrix6{i}=[data6{i}];matrix7{i}=[data7{i}];matrix8{i}=[data8{i}];matrix9{i}=[data9{i}];
end
d1=cell2mat(matrix1);d2=cell2mat(matrix2);d3=cell2mat(matrix3);d4=cell2mat(matrix4);d5=cell2mat(matrix5);
d6=cell2mat(matrix6);d7=cell2mat(matrix7);d8=cell2mat(matrix8);d9=cell2mat(matrix9);
for i=1:n
    temp1(i,:)=d1(((i-1)*m+1):(i*m));temp2(i,:)=d2(((i-1)*m+1):(i*m));temp3(i,:)=d3(((i-1)*m+1):(i*m));temp4(i,:)=d4(((i-1)*m+1):(i*m));
    temp5(i,:)=d5(((i-1)*m+1):(i*m));temp6(i,:)=d6(((i-1)*m+1):(i*m));temp7(i,:)=d7(((i-1)*m+1):(i*m));temp8(i,:)=d8(((i-1)*m+1):(i*m));
    temp9(i,:)=d9(((i-1)*m+1):(i*m));
end
temp_1=mean(temp1)/N;temp_4=mean(temp4)/N;temp_5=mean(temp5)/N;temp_6=mean(temp6)/N;temp_7=mean(temp7)/N;temp_8=mean(temp8)/N;
temp_9=mean(temp9)/N;temp_2=mean(temp2)/N;temp_3=mean(temp3)/N;%temp_3=mean(temp3)/N;
% temp_11=mean(temp11)/N;temp_44=mean(temp44)/N;temp_55=mean(temp55)/N;temp_66=mean(temp66)/N;temp_77=mean(temp77)/N;
% temp_88=mean(temp88)/N;temp_99=mean(temp99)/N;temp_1010=mean(temp1010)/N;%temp_22=mean(temp22)/N;temp_33=mean(temp33)/N;
%% 
figure();
space=time/16;
plot(temp_1(1:space:time),'ko-');hold on;plot(temp_2(1:space:time),'mp-');hold on;plot(temp_3(1:space:time),'r*-');hold on;
plot(temp_4(1:space:time),'gx-');hold on;plot(temp_5(1:space:time),'bs-');hold on;plot(temp_6(1:space:time),'cd-');hold on;
plot(temp_7(1:space:time),'^-','color',[1 0.8 0]);hold on;plot(temp_8(1:space:time),'<-','color',[0.2 0.5 1]);hold on;
plot(temp_9(1:space:time),'>-','color',[1 0.38 0]);hold on;
xlabel('infected time/t');ylabel('S(t)');
hleg2=legend('DC','KS','Cnc','Cnc+','C(ks)','\theta','MDD','MCDE','IKS');
set(hleg2, 'Location', 'SouthEast');
title('USAir');
%% 
figure();
step=floor(nlen/20);
x=linspace(0.01,0.2,20);
plot(x,I_mean_xianshi,'ko-');hold on;
plot(x,I2_mean_xianshi,'mp-');hold on;%(round([1:step:nlen]/step)+1)
plot(x,I3_mean_xianshi,'r*-');hold on;plot(x,I4_mean_xianshi,'gx-');hold on;
plot(x,I5_mean_xianshi,'bs-');hold on;plot(x,I6_mean_xianshi,'cd-');hold on;
plot(x,I7_mean_xianshi,'^-','color',[1 0.8 0]);hold on;
plot(x,I8_mean_xianshi,'<-','color',[0.2 0.5 1]);hold on;
plot(x,I9_mean_xianshi,'>-','color',[1 0.38 0]);hold on;
% axis([0,22,0,0.45]);
% set(gca,'XTickLabel',[0:0.02:0.2]);
xlabel('Proportion of initial infected nodes');ylabel('S');
hleg2=legend('DC','KS','Cnc','Cnc+','C(ks)','\theta','MDD','MCDE','IKS');
set(hleg2, 'Location', 'SouthEast');
title('USAir');
%% 
spreadtime=20;
for i=1:n
    spreading(i,:)=spreadingability(a,beta,mu,spreadtime);      
end
spreading1=mean(spreading);

figure();
subplot(1,3,1)
coeff1 = corr(DC', spreading1' , 'type' , 'Kendall') ;
scatter(DC',spreading1',10,'filled'); %filled
xlabel('DC');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff1))]);

coeff2 = corr(ks', spreading1' , 'type' , 'Kendall');
subplot(1,3,2);
scatter(ks',spreading1',10,'filled'); %filled
xlabel('KS');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff2))]);

coeff3 = corr(Cnc', spreading1' , 'type' , 'Kendall'); 
subplot(1,3,3);
scatter(Cnc',spreading1',10,'filled'); 
xlabel('Cnc');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff3))]);

figure();
subplot(1,3,1);
coeff4 = corr(Cncopt', spreading1' , 'type' , 'Kendall'); 
scatter(Cncopt',spreading1',10,'filled'); 
xlabel('Cnc+');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff4))]);

coeff5 = corr(C_fai', spreading1' , 'type' , 'Kendall'); 
subplot(1,3,2);
scatter(C_fai',spreading1',10,'filled'); 
xlabel('C(ks)');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff5))]);

coeff6 = corr(ksk', spreading1' , 'type' , 'Kendall'); 
subplot(1,3,3);
scatter(ksk',spreading1',10,'filled'); 
xlabel('\theta');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff6))]);

figure();
subplot(1,3,1);
coeff7 = corr(result_km', spreading1' , 'type' , 'Kendall'); 
scatter(result_km',spreading1',10,'filled'); 
xlabel('MDD');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff7))]);

coeff8 = corr(mcde', spreading1' , 'type' , 'Kendall'); 
subplot(1,3,2);
scatter(mcde',spreading1',10,'filled');
xlabel('ksh');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff7))]);

coeff9 = corr(e', spreading1' , 'type' , 'Kendall'); 
subplot(1,3,3);
scatter(e',spreading1',10,'filled'); 
xlabel('IKS');
ylabel('G_{i}');
title(['\tau = ',num2str(abs(coeff9))]);
L1=[];L2=[];L3=[];L4=[];L5=[];L6=[];L7=[];L8=[];L9=[];
for i=1:floor(nlen/20):nlen  %for i=1:62:nlen
   L1=[L1, sum(sum(d(z1(1:i),z1(1:i))))/(i*(i-1))]; 
   L2=[L2, sum(sum(d(z2(1:i),z2(1:i))))/(i*(i-1))];
   L3=[L3, sum(sum(d(zCnc(1:i),zCnc(1:i))))/(i*(i-1))];
   L4=[L4, sum(sum(d(zCncopt(1:i),zCncopt(1:i))))/(i*(i-1))];
   L5=[L5, sum(sum(d(zC_fai(1:i),zC_fai(1:i))))/(i*(i-1))];
   L6=[L6, sum(sum(d(zksk(1:i),zksk(1:i))))/(i*(i-1))];
   L7=[L7, sum(sum(d(zkm(1:i),zkm(1:i))))/(i*(i-1))];   
   L8=[L8, sum(sum(d(zmcde(1:i),zmcde(1:i))))/(i*(i-1))];
   L9=[L9, sum(sum(d(index(1:i),index(1:i))))/(i*(i-1))];   
end
figure()
temp=linspace(0.01,0.2,19)
plot(temp,L1(2:end),'ko-');hold on;
plot(temp,L2(2:end),'mp-');hold on;plot(temp,L3(2:end),'r*-');hold on;
plot(temp,L4(2:end),'gx-');hold on;plot(temp,L5(2:end),'bs-');hold on;plot(temp,L6(2:end),'cd-');hold on;
plot(temp,L7(2:end),'^-','color',[1 0.8 0]);hold on;plot(temp,L8(2:end),'<-','color',[0.2 0.5 1]);hold on;
plot(temp,L9(2:end),'>-','color',[1 0.38 0]);hold on;
axis([0,0.2,0,2.5]);
xlabel('Proportion of initial infected nodes');ylabel('L_s');
hleg2=legend('DC','KS','Cnc','Cnc+','C(ks)','\theta','MDD','MCDE','IKS');
set(hleg2, 'Location', 'southeast');title('USAir');

m_DC = evaluate_M(DC,N);
m_ks = evaluate_M(ks,N);
m_Cnc = evaluate_M(Cnc,N);
m_Cncopt = evaluate_M(Cncopt,N);
m_C_fai = evaluate_M(C_fai,N);
m_ksk = evaluate_M(ksk,N);
M_mdd;
m_ksh = evaluate_M(ksh,N);
m_mcde = evaluate_M(mcde,N);
m_myself = evaluate_M(result9,N);

toc