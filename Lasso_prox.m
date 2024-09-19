function [xp] = Lasso_prox(x0,eps,A,b,gamma)
%这是使用邻近梯度法求LASSO问题的函数 Lasso函数存储着函数信息的计算
% x0初始向量，256维，eps为精度
[fk,gk]=Lasso(x0,A,b,gamma);
xk=x0;
iter=0;
namb=eigs(A'*A,1);
ak=1/namb;
normg=1;
while abs(normg)>eps && iter<50000
    iter=iter+1;
    f1=fk;%记录当前点函数值
    %第一步先进行梯度下降
   dk=-gk;
   %计算下降方向   
   xk=xk+ak*dk; 
   max=zeros(256,1);

%进行完梯度下降线搜索，用临近算子做收缩
   si=sign(xk);
   for i=1:256
       m=abs(xk(i))-ak*gamma;
       if m>0
           xk(i)=si(i)*m;
       else
           xk(i)=0;
       end
   end
   [fk,gk]=Lasso(xk,A,b,gamma);
   normg=fk-f1;
    fprintf('在第%d次迭代，当前点函数值与上一点函数值差量为  %f\n',iter,normg);
end
xp=xk;
end