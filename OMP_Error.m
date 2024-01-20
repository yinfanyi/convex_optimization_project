function [A]=OMP_Error(D,X,Error_Goal); 
%输入：D：字典
%      X：表示信号
%      Error_Goal：最大允许误差
%输出：A：稀疏系数矩阵

[~,P]=size(X);
[n,~]=size(D);
E_2 = Error_Goal^2*n;
Max_Coefficient = n/2;
A = sparse(size(D,2),size(X,2));
for k=1:1:P
    x=X(:,k);
    Residual=x;
	Indx = [];
	Current_Residual_Norm_2 = sum(Residual.^2);
	j = 0;

    while Current_Residual_Norm_2>E_2 && j < Max_Coefficient
		j = j+1;
        Proj=D'*Residual;
        Pos=find(abs(Proj)==max(abs(Proj)));
        Pos=Pos(1);
        Indx(j)=Pos;
        a=pinv(D(:,Indx(1:j)))*x;
        Residual=x-D(:,Indx(1:j))*a;
		Current_Residual_Norm_2 = sum(Residual.^2);
   end
   if (~isempty(Indx))
       A(Indx,k)=a;
   end
end
return;
