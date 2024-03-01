function y = AdamsBashforth4(f,x,y0)
%
% AdamsBashforth4 uses the fourth-order Adams-Bashforth formula to solve
% a first-order initial-value problem in the form y' = f(x,y), y(x0) = y0.
%
% y = AdamsBashforth4(f,x,y0), where
%
%   f is an anonymous function representing f(x,y),
%   x is a vector representing the mesh points,
%   y0 is a scalar representing the initial value of y,
%
%   y is the vector of solution estimates at the mesh points.
%
y(1:4) = RK4(f,x(1:4),y0); n = length(x);
for i = 4:n-1,
    h = x(i+1)-x(i);
    y(i+1) = y(i)+h*(55*f(x(i),y(i))-59*f(x(i-1),y(i-1))+37*f(x(i-2),
y(i-2))-9*f(x(i-3),y(i-3)))/24;
end