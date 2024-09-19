function [xp] = FRmethod(x0,eps)
%x0是初始点,它应该是一个列向量，eps应取一个很小的数
%xp为极值点
[fk,gk,H]=exfun(x0);
normg=norm(gk);
xk=x0;
iter=0;
while normg>eps && iter<10000
    iter=iter+1;
    
    if iter==1
        dk=-gk;%定义初始的dk
    else
        bk=-dk'*H*(-gk)/(dk'*H*dk);%直接算法
        bk=normg^2/normg0^2;%通过残量计算
        dk=-gk+bk*dk;%计算新的下降方向
    end    
    
    ak=-gk'*dk/(dk'*H*dk);%计算步长
    xk=xk+ak*dk;%得到新的迭代点
    [fk,gk]=exfun(xk);
    normg0=normg;%记录旧的rk'*rk
    normg=norm(gk);
    fprintf('使用共轭梯度方法，在第%d次迭代，残量的范数为  %f\n',iter,normg);
    pause(0.5)
end
xp=xk;
end

