function [xp] = Newton(x0,eps)
%x0是初始点,它应该是一个列向量，eps应取一个很小的数
%xp为极值点
[fk,gk,H]=exfun(x0);
normg=norm(gk);
xk=x0;
iter=0;
while normg>eps && iter<10000
    iter=iter+1;
    dk=-inv(H)*gk;%计算下降方向
    ak=1;%步长
    xk=xk+ak*dk;
    [fk,gk]=exfun(xk);
    normg=norm(gk);
    fprintf('使用牛顿方法，在第%d次迭代，残量的范数为  %f\n',iter,normg);
    
end
xp=xk;
end

