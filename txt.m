function  [A,Du,txt_line_number] = txt(  )
%-- 将txt文件中的数据集改为--邻接矩阵--形式进行保存
%-- 转换为无向、无权图
%-- a为节点个数

fid = fopen('1USAir.txt','r');
txt_line_number=0;  %--txt文件中行的个数，每读一行，+1
a=332;%1133;%每个txt文件改这里
A=zeros(a,a); %--初始化矩阵，行列为节点个数%%邻接矩阵
Du = zeros(1,a);  %--初始化度值=0%%度值
T_line=zeros(1,3);   %--存储txt中每行分割后的3个元素值

while ~feof(fid)        % 判断是否为文件末尾   
    str_line=fgetl(fid);  %--读取txt文件中的一行 
    [token_1, remain_1] = strtok(str_line); %--取第一个元素的值        
    i=str2double(token_1); %---矩阵A对应行坐标的值   
    [token_2, remain_2] = strtok(remain_1);    
    j=str2double(token_2); %---矩阵A对应列坐标的值
    
    %****---矩阵A对应位置上的值赋1---****%
    A(i,j) = 1;
    A(j,i) = A(i,j); 
    %************************************% 
    
    Du(i) = sum(A(i,:));%--求节点的度值（因矩阵对称位置都赋值了，所以矩阵行、列的度值都要记录下来）
    Du(j) = sum(A(j,:));%--求节点的度值（否则，边的总和不是节点度值总和的2倍）
    
    txt_line_number=txt_line_number+1;
    
end   %---while

filename= sprintf('%s-a-%d.mat','shell',a);  %--- 存放到D盘的位置
save(filename,'a');

filename= sprintf('%s-a-%d.mat','shell_Du',a);
save(filename,'Du');

% sum(Deg);  %---正确时（sum(Deg)=2*txt_line_number）

disp('ok!');  

end