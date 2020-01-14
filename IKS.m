function [result_iks,index,result5,ze,e] = IKS(a,N,Du,ks)

for i=1:N
   p(i)=Du(i)/sum(Du);%1总的度2邻居节点的度
end
e=zeros(1,N); 
for i=1:N
    for j=1:N     
        if a(i,j)~=0      
        e(i)=e(i)-(1/ks(j))*p(j).*log(p(j));%这里用log而不是log2%e(i)=-sum(p(:,i).*log(p(:,i)));
        end
    end
%     e(i)=e(i)+1/(ks(i)+1)
end
[result5,ze] = sort(e,'descend');
index=[];
for j=1:400%几次迭代这个不重要
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
% index是选节点的顺序,然后值是按照节点熵来的