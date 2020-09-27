function alphak = line_search(xk, dk, func, grad)
alphak = 1;
c1 = 0.25;
c2 = 0.5;
eps = 0.3;
d1 = transpose(dk);
x1 = num2cell(xk + alphak*d1);
x2 = num2cell(xk + alphak*d1);
x3 = num2cell(xk);

while (func(x1{:}) > func(x3{:}) + (c1*alphak*transpose(grad(x3{:}))*dk)) && (abs(transpose(grad(x2{:}))*dk) > c2*(abs(transpose(grad(x3{:}))*dk)))
    alphak = alphak*eps;
    x1 = num2cell(xk + alphak*d1);
    x2 = num2cell(xk + alphak*d1);
end
end