%二次函数极小化问题
x=randn(5,1);
[xp1] = steepest(x,1e-5)%使用最速下降
[xp2] = Newton(x,1e-5)%使用牛顿法
[xp3] = FRmethod(x,1e-5)%使用共轭梯度方法
[xp4] = DFP(x,1e-5)%使用拟牛顿dfp方法
[xp5] = BFGS(x,1e-5)%使用拟牛顿bfgs方法
