function [xp,normg] = BFGSzengguang(xk,sigma,namb,eta0)
    [fk,gk] = Zengguang(xk,sigma,namb);
    normg=1;
    t=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %先使用一次Armoji线搜索
    dk=-gk;
    g1=gk;%保存好上一个gk
    x1=xk; %保存好上一个xk
    %使用Armijo线搜索确定ak
    ak=0.3;
    f0=fk;
    x0=xk;
    %先用固定步长代一遍
    xk=xk+ak*dk;
    [fk]=Zengguang(xk,sigma,namb);
    f1=fk;
    slope=dot(gk,dk);
    %验证Armoji条件，进行先搜索，ρ取0.5，c取0.3
    while f1>f0+0.3*ak*slope&&ak>0.01
      ak=ak/2;
       xk=x0+ak*dk;
       [fk,gk]=Zengguang(xk,sigma,namb);
       f1=fk;
    end
    normg=norm(fk);
   iter=0;
    %使用BFGS方法处理
    Hk=eye(2);
while normg>eta0 && iter<3000
    iter=iter+1;
    %计算拟牛顿方程的BFGS版Hk
    yk=gk-g1;
    v=xk-x1;
    Hk=(eye(2)-v*yk'./(v'*yk))*Hk*(eye(2)-yk*v'./(v'*yk))+v*v'./(v'*yk);
    %计算下降方向
    dk=-Hk*gk;
    %预留前一个gk，xk
    g1=gk;
    x1=xk;
    lamda1 = wolfe_powell(xk,dk,namb,sigma);
     
    xk=x1+lamda1*dk;
     
    %得到下一个xk下一行计算相关数据
    [fk,gk]=Zengguang(xk,sigma,namb);
    normg=norm(gk);
end
    xp=xk;    
end

