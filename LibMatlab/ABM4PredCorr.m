function y = ABM4PredCorr(f,x,y0)
%
% ABM4PredCorr uses the fourth-order Adams-Bashforth-Moulton predictor-
% corrector formula to solve a first-order initial-value problem in the
% form y' = f(x,y), y(x0) = y0.
%
% y = ABM4PredCorr(f,x,y0), where
%
% f is an anonymous function representing f(x,y),
% x is a vector representing the mesh points,
% y0 is a scalar representing the initial value of y,
%
% y is the vector of solution estimates at the mesh points.
%
py = zeros(4,1);
% Pre-allocate
y(1:4) = RK4(f,x(1:4),y0);
% Find the first 4 elements by RK4
h = x(2) - x(1); n = length(x);
% Start ABM4
for i = 4:n-1,
    py(i+1) = y(i) + (h/24)*(55*f(x(i),y(i))-59*f(x(i-1),y(i-1))+
37*f(x(i-2),y(i-2))-9*f(x(i-3),y(i-3)));
    y(i+1) = y(i) + (h/24)*(9*f(x(i+1),py(i+1))+19*f(x(i),y(i))-
5*f(x(i-1),y(i-1))+f(x(i-2),y(i-2)));
end