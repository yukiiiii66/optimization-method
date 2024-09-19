function [xp] = DFP(x0,eps)
%x0是初始点,它应该是一个列向量，eps应取一个很小的数
%xp为极值点
%步长搜索采用得是精确线方法
[fk,gk,H]=exfun(x0);
normg=norm(gk);
xk=x0;
iter=1;
Hk=eye(5);
dk=-Hk*gk;
ak=-gk'*dk/(dk'*H*dk);
x1=xk %保存好上一个xk
xk=xk+ak*dk;
g1=gk;%保存好上一个gk
[fk,gk,H]=exfun(xk);
normg=norm(gk);
fprintf('在第%d次迭代，残量的范数为  %f\n',iter,normg);
while normg>eps && iter<10000
    iter=iter+1;
    %计算拟牛顿方程的Hk
    yk=gk-g1;
    v=xk-x1;
    u=Hk*yk;
    alphak=-1/(yk'*Hk*yk);
    betak=1/(v'*yk);
    Hk=Hk+alphak*u*u'+betak*v*v';
    %计算下降方向
    dk=-Hk*gk;
    %计算步长
    ak=-gk'*dk/(dk'*H*dk);
     %预留前一个gk，xk
    g1=gk;
    x1=xk;
    %得到下一个xk
    xk=xk+ak*dk;
    [fk,gk]=exfun(xk);
    normg=norm(gk);
    fprintf('使用拟牛顿dfp算法，在第%d次迭代，残量的范数为  %f\n',iter,normg);
    pause(0.5)
end
xp=xk;
end

