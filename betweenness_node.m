function B=betweenness_node(A,a)
%%������ڵ������BY QiCheng
%%˼�룺�ڵ�i��j��ľ�����ڽڵ�i��k�������ڵ�k��j�����ʱ��i��j������·������k��
%%��Ϊi��j�ڵ�����·������kʱ��i��k��k��j�ض��������·����������÷�֤����֤����
% A�������������ڽӾ���������Ǹ�Ȩͼ
% a==0Ϊ�������磻a==1Ϊ�������磻
% B���������ڵ����
N=size(A,1);
B=zeros(1,N);
[D,C,aver_D]=Distance_F(A); %C��ij�����·������
if a==0
    for k=1:N
        for i=1:N
            if i~=k     %���˽ڵ㲻��               
                for j=i+1:N    %��Ϊ������ģ��������򡢷���ֻ��һ�Σ�����ֻ��һ�룻
                               %����Ļ��ۼ�һ����������ˣ�Ҳ��Ӱ��         
                    if j~=k    %���˽ڵ㲻��
                        if D(i,j)==D(i,k)+D(k,j)&C(i,j)~=0   %����������֤��ij�����·������k�ڵ�
                            B(k)=B(k)+C(i,k)*C(k,j)/C(i,j);
                        end
                    end
                end
            end
        end
    end
else
    for k=1:N
        for i=1:N
            if i~=k     %���˽ڵ㲻��
                for j=1:N
                    if j~=k      %���˽ڵ㲻��
                        if D(i,j)==D(i,k)+D(k,j)&C(i,j)~=0    %����������֤��ij�����·������k�ڵ�
                            B(k)=B(k)+C(i,k)*C(k,j)/C(i,j);
                        end
                    end
                end
            end
        end
    end
end
