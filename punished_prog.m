function [xp] = punished_prog(xk,eps,sigma)
%此为使用惩罚函数法处理不等式约束二次规划问题的程序
%需要调用punished_fun函数，它能给出二次罚函数的数据，注意初始点应给定应满足约束条件
sigma=sigma;
%ρ为罚因子增长系数
rou=1.2;iter=1;
[fk,gk] = punished_fun(xk,sigma);
%对初始数据进行一次存储
g1=gk;normg=1;
%对二次罚函数进行Armoji线搜索
    dk=-gk;%计算下降方向
   %计算步长
   ak=1;
   f0=fk;
   x0=xk;
   %先用固定步长代一遍
   xk=xk+ak*dk;
   [fk]=punished_fun(xk,sigma);
   f1=fk;
   slope=dot(gk,dk);
   %验证Armoji条件，进行先搜索，ρ取0.5，c取0.3
   while f1>f0+0.3*ak*slope
      ak=ak/2;
      xk=x0+ak*dk;
       [fk]=punished_fun(xk,sigma);
       f1=fk;
   end
   [~,gk]=punished_fun(xk,sigma);
   normg=norm(xk-x0);
    fprintf('在第%d次迭代，迭代点间距离为  %f\n',iter,normg);    
    %用增长系数更新罚因子
    sigma=rou*sigma;
%使用BB方法计算每一次的迭代步长完成后续循环
while normg>eps && iter<10000
    iter=iter+1;
    dk=-gk;%计算下降方向
   v=xk-x0;
    yk=gk-g1 ;   
    %使用BB方法计算步长
    ak=v'*yk./(yk'*yk)
    x0=xk;
    xk=xk+ak*dk;
    g1=gk;
    %计算新的点的数据
   [~,gk]=punished_fun(xk,sigma);
   normg=norm(xk-x0);
    fprintf('在第%d次迭代，两个迭代点的距离为  %f\n',iter,normg);
end

xp=xk;
end


