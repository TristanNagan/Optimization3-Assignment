function [xstar, output] = steepest_descent(f, g, x0, tol)
% === Inputs ===
% f   - function to be minimized
% g   - gradient of the function
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

%Initialise
xHist = x0;
iter = 0;
status = 1;

xcurr = x0;
xc = num2cell(xcurr);
fHist = f(xc{:});

while (norm(g(xc{:})) > tol) && (iter < 100000) % limiting the maximum number of iterations, no infinite loop 
    xcurr = xHist(end, :);
    xc = num2cell(xcurr);
    
    %Search Direction
    d = -g(xc{:});
    if (d == 0)
        break;
    end
    
    %calculating alpha by means of LS using strong wolfe's condition
    alpha = line_search2(xcurr, d, f, g);
    
    % updating x(k+1)
    xnew = xcurr + transpose(alpha*d);
    xn = num2cell(xnew);
    
    % updating the xHist as well as the fHist
    xHist = [xHist; xnew];
    fHist = [fHist; f(xn{:})];
    
    % updating the iteration count
    iter = iter + 1; 
end

% checking to see if the grad reached tolerence
if(norm(g(xc{:})) < tol) 
    status = 0;
end

xstar = xcurr;
xs = num2cell(xstar);
output = struct('iter', iter, 'status', status, 'fstar', f(xs{:}), 'xHist', xHist, 'fHist', fHist);
fprintf('xstar: [%d, %d]\n', xstar(1), xstar(2));
disp(output);
end