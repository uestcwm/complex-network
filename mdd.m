function [result_km,M,km,res_index,append]=mdd(a,NodeNum)%a���ڽӾ���
lamda=0.7;
%  a=[0 1 0 0 1 1 0 0;
%     1 0 0 0 1 1 0 0;
%     0 0 0 0 0 1 0 0;
%     0 0 0 0 0 0 1 0
%     1 1 0 0 0 1 0 0;
%     1 1 1 0 1 0 1 0;
%     0 0 0 1 0 1 0 1;
%     0 0 0 0 0 0 1 0];
%{

%}
% NodeNum=size(a,2);
du=sum(a);
km=du;
kr=du;
ki=zeros(1,NodeNum);
append=[];
res_index=[];
result_km=zeros(1,NodeNum);
while sum(sum(a))~=0
        [col]=find(km==min(km(km>0)));%row��col��
        res_index=[res_index;min(km(km>0))];
        append=[append;length(col)];
        result_km([col])=min(km(km>0));        
        for j=1:length(col)
            length(col);
            [row]=find(a(col(j),:)~=0);%��һ��3�����ӵ�
            ki(row)=ki(row)+1;    
        end
        km=zeros(1,NodeNum);%ÿһ��Ҫ��������
        a(col,:)=0;
        a(:,col)=0;%�൱�ڰ���������ɾ��
        ki(col)=0;
        kr=sum(a);%����ͼ�еĽڵ�Ķ�   
        km=kr+0.7*ki; 
end  
[col]=find(km==min(km(km>0)));
result_km([col])=min(km(km>0));
res_index=[res_index;min(km(km>0))];
app=length(find(km~=0));
append=[append;app];
sumer=0;
for i =1:length(append)
    sumer=sumer+append(i)*(append(i)-1);
end
M=(1-sumer/(NodeNum*(NodeNum-1)))^2;
end