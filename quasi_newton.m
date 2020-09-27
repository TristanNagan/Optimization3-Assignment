function [xstar, output] = quasi_newton(func, grad, x0, tol)
% === Inputs ===
% func   - function to be minimized
% grad   - gradient of the function
% x0  - initial solution
% tol - tolerance value(for stopping)

% === Outputs ===
% x        - optimal solution
% output.  - matlab structure with the following fields
%     .iter   - number of iteration
%     .status - 0: gradient reached tolerance
%               1: reached maximum number of iterations
%     .fstar  - optimal function value f(x*)
%     .xHist  - history of x i.e store x at each iteration.
%     .fHist  - history of f i.e store f at each iteration.

% func -> initial size of the Hessian for the different functions 
sz = size(grad);
H = eye(sz(1));

xHist = x0;
iter = 0;
status = 1;
xcurr = x0;
xc = num2cell(xcurr);
fHist = func(xc{:});
i = 0;
aHist = 1;


while (norm(grad(xc{:})) > tol)  && (iter < 100000)
    xcurr = xHist(end, :);
    xc = num2cell(xcurr);
    if ((grad(xc{:})) == 0)
        break;
    end
    qu = grad(xc{:});
    d = -H*grad(xc{:}); % setting the search direction, should a 2x1
    i = i + 1;
    if(isnan(d))
        disp('wtf');
    end
    if iter ~= 0
        alph = line_search2(xcurr, d, func, grad);
    else
        alph = 1;
    end
    %alph = 1;
    
    aHist = [aHist; alph];
    
    xnew = xcurr + transpose(alph*d);  %xcurr is row vector, alpha*d is column so cant add
    xn = num2cell(xnew);
    
    xHist = [xHist; xnew];
    fHist = [fHist; func(xn{:})];
        % compute g and update H

    gam = grad(xn{:})- grad(xc{:}); % column vector, gamma
    del = (xnew - xcurr);  % column vector, delta
    
    % breaking up the update into parts
    t1 = (transpose(del)*del)/(del*gam);
    t2 = (H*gam*transpose(H*gam))/(transpose(gam)*H*gam);
    z = (transpose(del))/(del*gam)-(H*gam)/(transpose(gam)*H*gam);
    t3 = transpose(gam)*H*gam*z*transpose(z);
    H = H + t1 + t2 + t3; % updating H
    iter = iter + 1; % updating iterations
end
if(norm(grad(xc{:})) < tol)
    status = 0;
end
xstar = xcurr;
xs = num2cell(xstar);
output = struct('iter', iter, 'status', status, 'fstar', func(xs{:}), 'xHist', xHist, 'fHist', fHist);
fprintf('xstar: [%d, %d]\n', xstar(1), xstar(2));
disp(output);
end