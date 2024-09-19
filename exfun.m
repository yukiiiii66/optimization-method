function [f,g,H] = exfun(x)
%x 为输入点，输出结果为函数值，梯度以及海森阵信息
%   x 是五维向量
f=0;
s=0;
for i=1:5
    for j=1:5
      s=s+x(j);
    end
   
    f=f+5*x(i)^2+i*x(i)+x(i)*s;
     s=0;
end
g=zeros(5,1);
for i=1:5
    g(i)=12*x(i)+2*(sum(x)-x(i))+i;
end
H=[12,2,2,2,2;2,12,2,2,2;2,2,12,2,2;2,2,2,12,2;2,2,2,2,12];
end
