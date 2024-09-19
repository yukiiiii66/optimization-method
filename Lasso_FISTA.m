function [xp] = Lasso_FISTA(x0,eps,A,b,gamma)
%这是使用加速邻近梯度法求LASSO问题的函数 Lasso函数存储着函数信息的计算
% x0初始向量，256维，eps为精度
[fk,gk]=Lasso(x0,A,b,gamma);
normg=1;
xk=x0;
iter=0;
namb=eigs(A'*A,1);
ak=1/(1.5*namb);
vk=xk;
while abs(normg)>eps && iter<50000
    iter=iter+1;
    f1=fk; 
    %计算加速后的点在这个点基础上做刚才的临近梯度方法
    gam=2/(iter+1);
    %使用改进的FISTA算法，gam这个系数会随着迭代次数二改变
    yk=(1-gam)*xk+gam*vk;
    x1=xk;
    %第一步先进行梯度下降
    dk=-gk;
   xk=yk+ak*dk; 
%进行完梯度下降线搜索后，用临近算子做收缩
   si=sign(xk);
   for i=1:256
       m=abs(xk(i))-ak*gamma;
       if m>0
           xk(i)=si(i)*m;
       else
           xk(i)=0;
       end
   end
   %计算下一次迭代需要的vk
   vk=x1+(1/gam)*(xk-x1);
   [fk,gk]=Lasso(xk,A,b,gamma);
   normg=fk-f1;
    fprintf('在第%d次迭代，当前点函数值与上一点函数值差量为  %f\n',iter,normg);
end
xp=xk;
end