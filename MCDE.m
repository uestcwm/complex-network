function mcd = MCDE( a,du,ks,N )
% du=sum(a);
% ks=mycoreness(a);
% N=size(a,1);
mcd=zeros(1,N);
entropy=zeros(1,N);
for i=1:N
  for j=1:max(ks)
    temp=sum(ks(find(a(i,:)~=0))==j);
    if temp==0
        continue
    end
    entropy(i)=-(temp/du(i))*(log(temp/du(i))/log(2))+entropy(i);
  end
  mcd(i)=entropy(i)+du(i)+ks(i);
end
end