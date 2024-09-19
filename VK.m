function [v] = VK(xk,sigma,namb)
%此为在增广拉格朗日函数中关于约束违反度的计算程序，返回约束违反度

v1=-11+xk(1)^2-6*xk(1)+4*xk(2);
if v1>-namb(1)/sigma
    v1=v1^2;
else
    v1=(-namb(1)/sigma)^2;
end

v2=-xk(1)*xk(2)+3*xk(2)+exp(xk(1)-3)-1;
if v2>-namb(2)/sigma
    v2=v2^2;
else
    v2=(-namb(2)/sigma)^2;
end
v3=-xk(1);
if v3>-namb(3)/sigma
    v3=v3^2;
else
    v3=(-namb(3)/sigma)^2;
end
v4=-xk(2);
if v4>-namb(4)/sigma
    v4=v4^2;
else
    v4=(-namb(4)/sigma)^2;
end
v=sqrt(v1+v2+v3+v4);

end

