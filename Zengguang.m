function [f,g] = Zengguang(x,sigma,namb)
%此为该二次规划问题增广拉格朗日函数，主要求其关于x的梯度，在给定乘子向量namb(4维向量)与x下的函数值

f0=x(1)^2+x(2)^2-16*x(1)-10*x(2);
f1=-11+x(1)^2-6*x(1)+4*x(2);
%计算第一个增广项
if namb(1)/sigma+f1>0
    f1=namb(1)/sigma+f1;
    g1=2*(namb(1)/sigma+f1).*[2*x(1)-6;4];
else
    f1=0;
    g1=[0;0];
end
f1=f1^2-(namb(1)^2)/(sigma^2);
%计算第二个增广项
f2=-x(1)*x(2)+3*x(2)+exp(x(1)-3)-1;
if namb(2)/sigma+f2>0
    f2=namb(2)/sigma+f2;
    g2=2*(namb(2)/sigma+f2).*[-x(2)+exp(x(1)-3);-x(1)+3];
else
    f2=0;
    g2=[0;0];
end
f2=f2^2;
%计算第三个增广项
f3=-x(1);
if namb(3)/sigma+f3>0
    f3=namb(3)/sigma+f3;
    g3=2*namb(3)/sigma+f3.*[-1;0];
else
    f3=0;
    g3=[0;0];
end
f3=f3^2-(namb(3)^2)/(sigma^2);
%计算第四个增广项
f4=-x(2);
if namb(4)/sigma+f4>0
    f4=namb(4)/sigma+f4;
    g4=2*namb(4)/sigma+f4.*[0;-1];
else
    f4=0;
    g4=[0;0];
end
f4=f4^2-(namb(4)^2)/(sigma^2);
f=f0+(sigma/2)*(f1+f2+f3+f4);
%计算关于x的梯度
g=[2*x(1)-16;2*x(2)-10]+(sigma/2).*(g1+g2+g3+g4);

end

