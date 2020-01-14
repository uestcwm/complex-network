function W = CChole( a,SH,CC )
du=sum(a);
N=size(a,2);
W=zeros(1,N); 
for i=1:N
    for j=1:N     
        if a(i,j)~=0      
        W(i)=W(i)+CC(j)/SH(j);%这里用log而不是log2%e(i)=-sum(p(:,i).*log(p(:,i)));
        end
    end
end
for i=1:N
    for j=1:N
            W(i)=(W(i)+1/SH(i))/N;%邻居度   
    end
end

end

