function [result_iks,index,result5,ze,e] = IKS(a,N,Du,ks)

for i=1:N
   p(i)=Du(i)/sum(Du);%1�ܵĶ�2�ھӽڵ�Ķ�
end
e=zeros(1,N); 
for i=1:N
    for j=1:N     
        if a(i,j)~=0      
        e(i)=e(i)-(1/ks(j))*p(j).*log(p(j));%������log������log2%e(i)=-sum(p(:,i).*log(p(:,i)));
        end
    end
%     e(i)=e(i)+1/(ks(i)+1)
end
[result5,ze] = sort(e,'descend');
index=[];
for j=1:400%���ε����������Ҫ
    k=max(ks);
    for k=k:-1:1
       [row,column]=find(ks==k);
       [row2,value]=max(e(column));
       index=[index;column(value)];
    end
    ks(index)=0;
end
result_iks=e(index);
end
% index��ѡ�ڵ��˳��,Ȼ��ֵ�ǰ��սڵ�������