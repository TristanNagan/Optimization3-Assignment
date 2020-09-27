function alphak = line_search2(xk, dk, func, grad)
alpham = 20;
alphap = 0;
c1 = 1e-4;
c2 = 0.9;
alphax = 1;
maxIter = 100;
i = 1;
%gx0 = gx0'*d;
%fxp = fx0;
%gxp = gx0;

while true
    d = transpose(dk);
    xx = num2cell(xk + alphax*d);
    x = num2cell(xk);
    if func(xx{:}) > func(x{:}) + c1*alphax*transpose(grad(x{:}))*dk
        alphak = zoom(func, grad, xk, dk, alphap, alphax);
        return;
    end
    if abs(grad(xx{:})) <= -c2*grad(x{:})*d
        alphak = alphax;
        return;
    end
    if grad(xx{:}) >= 0
        alphak = zoom(func, grad, xk, dk, alphax, alphap);
        return;
    end
    
    alphap = alphax;
    
    if i > maxIter
        alphak = alphax;
        return;
    end
    
    r = 0.8;
    alphax = alphax + (alpham - alphax)*r;
    i = i + 1;
end

end