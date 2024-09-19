n=256;m=128;
A=randn(m,n);u=sprandn(n,1,0.1);
b=A*u;gamma=10;
disp('使用邻近梯度法处理lasso问题')
[xp1] = Lasso_prox(ones(256,1),1e-8,A,b,gamma);
%使用Nesterov加速邻近梯度法处理lasso问题
[xp2] = Lasso_FISTA(ones(256,1),1e-8,A,b,gamma);
%使用交替乘子法处理Lasso问题
[xp3] = Lasso_ADMM(ones(256,1),1e-8,A,b,gamma);

normg=norm(xp1-xp2)
normg=norm(xp2-xp3)
normg=norm(xp3-xp1)