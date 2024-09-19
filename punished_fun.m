function [f,g] = punished_fun(x,sigma)
%此为不等式约束最优化问题的二次罚函数
%   sigama为罚因子
g1=-1*(11-x(1)^2+6*x(1)-4*x(2));
if g1>0
   c1=g1^2;
   grad1=[(12-4*x(1))*(11-x(1)^2+6*x(1)-4*x(2));-8*(11-x(1)^2+6*x(1)-4*x(2))];
elseif g1<=0
    c1=0;
    grad1=[0;0];
end
g2=-1*(x(1)*x(2)-3*x(2)-exp(x(1)-3)+1);
if g2>0
    c2=g2^2;
    grad2=[2*(x(1)*x(2)-3*x(2)-exp(x(1)-3)+1)*(x(2)-exp(x(1)-3));2*(x(1)*x(2)-3*x(2)-exp(x(1)-3)+1)*(x(1)-3)];
elseif g2<=0
    c2=0;
    grad2=[0;0];
end
g3=-x(1);
if g3>0
    c3=g3^2;
    grad3=[2*g3;0];
elseif g3<=0
    c3=0;
    grad3=[0;0];
end
g4=-x(2);
if g4>0
    c4=g4^2;
    grad4=[0;2*g4];
elseif g4<=0
    c4=0;
    grad4=[0;0];
end
f=x(1)^2+x(2)^2-16*x(1)-10*x(2)+0.5*sigma*(c1+c2+c3+c4);
g=[2*x(1)-16;2*x(2)-10]+0.5*sigma.*(grad1+grad2+grad3+grad4);
end

