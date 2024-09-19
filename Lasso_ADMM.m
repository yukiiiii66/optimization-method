function [xp] = Lasso_ADMM(x0,eps,A,b,gamma)
%此算法为交替乘子算法，对于该方法我们是将无约束优化问题转换为一个带线性约束可分离的优化问题
%具体采用的是最优化建模教材P473页介绍的Douglas-Rachford Splitting算法
[fk,~]=Lasso(x0,A,b,gamma);
%xk，zk是原问题所转化的分离优化问题的两个变量，满足Ax-b=z,f是Lasso问题函数的第一项，h是Lasso问题函数第二项
xk=x0;
iter=0;
normg=1;
t=0.5;
zk=xk;
while (normg>eps||normg<=0)&&iter<10000
%使用t*h的邻近算子计算xk,t是正常数
iter=iter+1;
f1=fk;
si=sign(zk);
   for i=1:256
       m=abs(zk(i))-t*gamma;
       if m>0
           xk(i)=si(i)*m;
       else
           xk(i)=0;
       end
   end
%计算yk 由光滑部分tf的邻近算子计算
yk=2*xk-zk;
yk=inv(t*A'*A+eye(256))*(t*A'*b+yk);
%交替乘子思想再计算下一个zk
zk=zk+yk-xk;
%最后计算两次迭代函数间的差值
[fk,~]=Lasso(xk,A,b,gamma);
normg=fk-f1;
fprintf('在第%d次迭代，当前点函数值与上一点函数值差量为  %f\n',iter,normg);
end
xp=xk;
end

