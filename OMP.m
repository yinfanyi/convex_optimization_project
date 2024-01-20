function [A]=OMP(D,X,L)
%输入： 
%    D：列归一化的字典
%    X：表征信号
%    L：各信号最大系数值
%输出： 
%    A：稀疏系数矩阵

[~,P]=size(X);
[~,K]=size(D);
for k=1:1:P
    a=[];
    x=X(:,k);
    Residual=x;
    Indx=zeros(L,1);
    for j=1:1:L
        Proj=D'*Residual;
        [~,pos]=max(abs(Proj));
        pos=pos(1);
        Indx(j)=pos;
        a=pinv(D(:,Indx(1:j)))*x;
        Residual=x-D(:,Indx(1:j))*a;
        if sum(Residual.^2) < 1e-6
            break;
        end
    end
    Temp=zeros(K,1);
    Temp(Indx(1:j))=a;
    A(:,k)=sparse(Temp);
end
return;
