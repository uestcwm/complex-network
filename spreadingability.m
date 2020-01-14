function spread=spreadingability(A,beta,mu,time)
for i=1:length(A)
    InitialStatel=zeros(length(A),1);
    InitialStatel(i)=1;
    [I_spread,R_spread]=sire(A,InitialStatel,beta,mu,time);
    spread(i)=(I_spread(time)+R_spread(time))/length(A);%节点i的传播能力
end
end

