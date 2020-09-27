function H = Hess(f)
s = symvar(f);
H = hessian(f, s);
end