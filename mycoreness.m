function core = mycoreness( a )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
du=sum(a);
dminmax = minmax(du);
core=dminmax(1)*ones(size(du));
td=du;
for k=dminmax(1):dminmax(2)
    while sum(td>=1&td<=k)
        ind=find(td>=1&td<=k);
        a(:,ind)=0;a(ind,:)=0;
        td=sum(a);
    end
    core=core+(td>0);
end

