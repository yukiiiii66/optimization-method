function lamda1 = wolfe_powell(xk,dk,namb,sigma)
c1 = 0.1;c2=0.5;
a = 0; b =Inf;
lamda1 = 1;
while(1)
    lamda1;
    [~, g] = Zengguang(xk,sigma,namb);
    if lamda1<0.001
        break;
    end
    if ~(Zengguang(xk+lamda1*dk,sigma,namb)-Zengguang(xk,sigma,namb) <= c1*g'*dk)
        b = lamda1;
        lamda1 = (lamda1 + a)/2;
        continue;
    end
    [~, g2] = Zengguang(xk+lamda1*dk,sigma,namb);
    if ~((g2)'*dk >= c2*(g)'*dk)
        a = lamda1;
        lamda1 = min([2*lamda1,(b+lamda1)/2]);
        continue;
    end
    break;
end
end

