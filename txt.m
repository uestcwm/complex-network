function  [A,Du,txt_line_number] = txt(  )
%-- ��txt�ļ��е����ݼ���Ϊ--�ڽӾ���--��ʽ���б���
%-- ת��Ϊ������Ȩͼ
%-- aΪ�ڵ����

fid = fopen('1USAir.txt','r');
txt_line_number=0;  %--txt�ļ����еĸ�����ÿ��һ�У�+1
a=332;%1133;%ÿ��txt�ļ�������
A=zeros(a,a); %--��ʼ����������Ϊ�ڵ����%%�ڽӾ���
Du = zeros(1,a);  %--��ʼ����ֵ=0%%��ֵ
T_line=zeros(1,3);   %--�洢txt��ÿ�зָ���3��Ԫ��ֵ

while ~feof(fid)        % �ж��Ƿ�Ϊ�ļ�ĩβ   
    str_line=fgetl(fid);  %--��ȡtxt�ļ��е�һ�� 
    [token_1, remain_1] = strtok(str_line); %--ȡ��һ��Ԫ�ص�ֵ        
    i=str2double(token_1); %---����A��Ӧ�������ֵ   
    [token_2, remain_2] = strtok(remain_1);    
    j=str2double(token_2); %---����A��Ӧ�������ֵ
    
    %****---����A��Ӧλ���ϵ�ֵ��1---****%
    A(i,j) = 1;
    A(j,i) = A(i,j); 
    %************************************% 
    
    Du(i) = sum(A(i,:));%--��ڵ�Ķ�ֵ�������Գ�λ�ö���ֵ�ˣ����Ծ����С��еĶ�ֵ��Ҫ��¼������
    Du(j) = sum(A(j,:));%--��ڵ�Ķ�ֵ�����򣬱ߵ��ܺͲ��ǽڵ��ֵ�ܺ͵�2����
    
    txt_line_number=txt_line_number+1;
    
end   %---while

filename= sprintf('%s-a-%d.mat','shell',a);  %--- ��ŵ�D�̵�λ��
save(filename,'a');

filename= sprintf('%s-a-%d.mat','shell_Du',a);
save(filename,'Du');

% sum(Deg);  %---��ȷʱ��sum(Deg)=2*txt_line_number��

disp('ok!');  

end