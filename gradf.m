function g = gradf(f)
s = symvar(f);
g = gradient(f, s);
end