function theta=ks_k(a,N,d,ks) %a���ڽӾ��� N�ǽڵ���
% [a,Du,txt_line_number] = txt(  );
%  a=[0 1 0 1 1 0;
%     1 0 0 1 1 0;
%     0 0 0 0 1 0;
%     1 1 0 0 1 0;
%     1 1 1 1 0 1;
%     0 0 0 0 1 0];
% N=size(a,2);
% BB=sparse(a);%ϡ�账��
% d=graphallshortestpaths(BB,'Directed',0);
% ks=mycoreness(a);%ascend
ks_max=max(ks);
[row,col]=find(ks==max(ks));%�ҵ����ksֵ����
for i=1:N
    distance(i)=sum(d(i,col));
    theta(i)=(ks_max-ks(i)+1)*distance(i);
end

