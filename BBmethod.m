function [xp] = BBmethod(x0,eps)
%BB方法，我们选择BB1格式
[fk,gk] = FUN2information(x0);
xk=x0;iter=0;
x1=xk;g1=gk;dk=-gk;
%使用Armijo线搜索确定x1以便后续计算
ak=1;
f0=fk;x0=xk;
%先用固定步长代一遍
   xk=xk+ak*dk;
   [fk]=FUN2information(xk);
   f1=fk;
   slope=dot(gk,dk);
 %验证Armoji条件，进行先搜索，ρ取0.5，c取0.3
   while f1>f0+0.3*ak*slope
      ak=ak/2;
       xk=x0+ak*dk;
       [fk]=FUN2information(xk);
       f1=fk;
   end
   [~,gk]=FUN2information(xk);
normg=norm(gk);
fprintf('在第%d次迭代，残量的范数为  %f\n',iter,normg);
%后面的ak采用BB方法进行计算
while normg>eps && iter<10000
    iter=iter+1;
    dk=-gk;%计算下降方向
    %计算BB方法下的步长ak
    yk=gk-g1;
    s=xk-x1;
    ak1=s'*s/(s'*yk);
    ak2=s'*yk/(yk'*yk);
    if ak1<ak2
        ak=ak1;
    else
        ak=ak2;
    end
    x1=xk;
    g1=gk;
    xk=xk+ak*dk;
    [~,gk]=FUN2information(xk);
   normg=norm(gk);
    fprintf('在第%d次迭代，残量的范数为  %f\n',iter,normg);
end
xp=xk;
end
