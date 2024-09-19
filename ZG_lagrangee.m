function [xp] = ZG_lagrangee(x0,sigma,namb,eps,eta,a,b,rou)
%此为增广拉格朗日法求解二次规划的函数程序 其中x0为初始点，必须满足约束条件，sigma为惩罚因子最好大于1
%namb为针对四个约束条件的四维乘子向量，eps为约束违反常数，eta为精度常数
%0<a<=b<1，rou>1为罚因子的增长系数
eta0=1/sigma;eps0=1/(sigma)^a;
%以给定的的x0为初始点进行Armoji线搜索求下一个迭代点
xk=x0;
[fk,gk] = Zengguang(xk,sigma,namb);
for iter=1:30
   
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
  while normg>eta0 && iter<13000
     iter=iter+1;
     %计算拟牛顿方程的BFGS版Hk
     yk=gk-g1;
     v=xk-x1;
     Hk=(eye(2)-v*yk'./(v'*yk))*Hk*(eye(2)-yk*v'./(v'*yk))+v*v'./(v'*yk);
     %计算下降方向
     dk=-Hk*gk
     %预留前一个gk，xk
     g1=gk;
     x1=xk;
     %步长由wolfe-powell准则确定
     lamda1 = wolfe_powell(xk,dk,namb,sigma);
     
     xk=x1+lamda1*dk;
     
     %得到下一个xk下一行计算相关数据
     [fk,gk]=Zengguang(xk,sigma,namb);
     normg=norm(gk);
  end
xk
  %计算约束违反度
  [v]=VK(xk,sigma,namb);
 if v<=eps0
      if v<=eps&&normg<=eta
          break;
      else
      end
      %惩罚因子不变，更新乘子
      s1=namb(1)+sigma*(-11+xk(1)^2-6*xk(1)+4*xk(2));
      if s1>0
          namb(1)=s1;
      else
          namb(1)=0;
      end
      s2=namb(2)+sigma*(-xk(1)*xk(2)+3*xk(2)+exp(xk(1)-3)-1);
      if s2>0
          namb(2)=s2;
      else
          namb(2)=0;
      end
      s3=namb(3)+sigma*(-xk(1));
      if s3>0
          namb(3)=s3;
      else
          namb(3)=0;
      end
      s4=namb(4)+sigma*(-xk(2));
       if s4>0
          namb(4)=s4;
      else
          namb(4)=0;
       end 
       %乘子更新完，更新一遍求解误差与约束违反度
      eta0=eta0/sigma;
      eps0=eps0/(sigma)^b;
 else
      %乘子不变更新，更新惩罚因子
      sigma=rou*sigma;
      %重新设置求解误差与约束违反度
      eta0=1/sigma;
      eps0=1/(sigma)^a;
 end
xp=xk;
end
end

