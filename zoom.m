function alphas = zoom(f, g, xk, dk, alphal, alphah)
maxIter = 100;
i = 0;
c1 = 1e-4;
c2 = 0.9;

while true
    alphax = 0.5*(alphal+alphah);
    alphas = alphax;
    d = transpose(dk);
    xx = num2cell(xk + alphax*d);
    xl = num2cell(xk + alphal*d);
    x = num2cell(xk);
    if (f(xx{:}) > f(x{:}) + c1*alphax*g(x{:})*d) | (f(xx{:}) >= f(xl{:}))
        alphah = alphax;
    else
        if abs(g(xx{:})) <= -c2*g(x{:})*d
            alphas = alphax;
            return;
        end
        if g(xx{:})*(alphah - alphal) >= 0
            alphah = alphal;
        end
        alphal = alphax;
    end
    i = i+1;
    if i > maxIter
        alphas = alphax;
        return;
    end

end
end