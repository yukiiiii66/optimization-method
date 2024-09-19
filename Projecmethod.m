function [xp] = Projecmethod(x0,eps)
%本题的约束区域是一个闭矩形区域，请确保你的起始点在[-1,0]x[1,2]内
[fk,gk] = FUN2information(x0);
normg=norm(gk);
xk=x0;
iter=0;
while abs(normg)>eps&&iter<200000
    iter=iter+1;
    %计算下降方向，做一次Armoji线搜索的梯度下降法
    dk=-gk;
   %计算步长
   ak=1;
   f0=fk;
   x0=xk;
   %先用固定步长代一遍
   xk=xk+ak*dk;
   [fk]=FUN2information(xk);
   f1=fk;
   slope=dot(gk,dk);
   %验证Armoji条件，进行先搜索，ρ取0.333，c取0.3
   while f1>f0+0.3*ak*slope
      ak=ak/3;
       xk=x0+ak*dk;
       [fk]=FUN2information(xk);
       f1=fk;
   end
   %梯度下降后，计算投影，来求得新的迭代点
   a1=[-1;1];b1=[0;2];
   for i=1:2
          if xk(i)<=a1(i)
               xk(i)=a1(i);
          elseif xk(i)>=b1(i)
               xk(i)=b1(i);
          end
   end       
   [fk,gk]=FUN2information(xk);
   normg=fk-f0;
    fprintf('在第%d次迭代，当前点函数值与上一点函数值差量为  %f\n',iter,normg);
end
xp=xk;
end


