function [f,g] = Lasso(x,A,b,gamma)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
f=1/2*norm(A*x-b)+gamma*norm(x,1);
g=A'*(A*x-b);
end