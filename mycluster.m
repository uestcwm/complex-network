function [ c ] = mycluster( a )
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
n=length(a);
for i=1:n
    m=find(a(i,:));
    ta=a(m,m);
    Lta=tril(ta);
    if length(m)==0|length(m)==1
        c(i)=0;
    else
        c(i)=sum(sum(Lta))/nchoosek(length(m),2);
    end
end
end

