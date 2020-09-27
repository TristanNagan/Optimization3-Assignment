function test()
syms x1 x2 x3 x4

%Functions:
Rosenbrock = 100*(x2-x1^2)^2+(1-x1)^2;
Powell = (x1+10*x2)^2+5*(x3-x4)^2+(x2-2*x3)^4+10*(x1-x4)^4;
Matyas = 0.26*(x1^2+x2^2)-0.48*x1*x2;

%Initial Points
initialPointsR = {[-1, 1]; [-3, -4]; [-2, 2]};%For Rosenbrock (Solution is [1, 1])
initialPointsP = {[5, 5, 5, 5]; [-4, -2, 2, 4]; [-1, 0, -3, 5]};%For Powell (Solution is [0, 0, 0, 0])
initialPointsM = {[-10, 10], [1, -1], [4, 5]};%For Matyas (Solution is [0, 0])

%Methods Performed on Rosenbrock
Grx = gradf(Rosenbrock);
Fxr = matlabFunction(Rosenbrock);
Gxr = matlabFunction(Grx);
Hr = Hess(Rosenbrock);
disp('Rosebrock Function')
for i = 1:3
    xri = initialPointsR(i,:);
    xr = cell2mat(xri);
    fprintf('Initial point: [%d, %d]\n', xr(1), xr(2));
    disp('Steepest Descent');
    SDr = steepest_descent(Fxr, Gxr, xr, 0.000001);
%   NMr = newton(Fxr, Gxr, Hr, xr, tol);
%   CGMr = cga(Fxr, Gxr, xr, tol);
    disp('QNM');
    [QNMr, out] = quasi_newton(Fxr, Gxr, xr, 0.000001);
end

%Methods Performed on Powell
Gpx = gradf(Powell);
Fxp = matlabFunction(Powell);
Gxp = matlabFunction(Gpx);
Hp = Hess(Powell);
disp('Powell Function');
for i = 1:3
    xpi = initialPointsP(i,:);
    xp = cell2mat(xpi);
    fprintf('Initial point: [%d, %d, %d, %d]\n', xp(1), xp(2), xp(3), xp(4));
    disp('Steepest Descent');
    SDp = steepest_descent(Fxp, Gxp, xp, 0.000001);
%   NMp = newton(Fxp, Gxp, Hp, xp, tol);
%   CGMp = cga(Fxp, Gxp, xp, tol);
    disp('QNM');
    QNMp = quasi_newton(Fxp, Gxp, xp, 0.000001);
end

%Methods Performed on Matyas
Gmx = gradf(Matyas);
Fxm = matlabFunction(Matyas);
Gxm = matlabFunction(Gmx);
Hm = Hess(Matyas);
disp('Matyas Function');
for i = 1:1
    xmi = initialPointsM(i,:);
    xm = cell2mat(xmi);
    fprintf('Initial point: [%d, %d]\n', xm(1), xm(2));
    disp('Steepest Descent');
    SDm = steepest_descent(Fxm, Gxm, xm, 0.000001);
%   NMm = newton(Fxm, Gxm, Hm, xm, tol);
%   CGMm = cga(Fxm, Gxm, xm, tol);
    disp('QNM');
    QNMm = quasi_newton(Fxm, Gxm, xm, 0.000001);
end

end