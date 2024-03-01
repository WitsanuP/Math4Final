function y = EulerODE(f,x,y0)
%
% EulerODE uses Euler's method to solve a first-order initial-value
% problem in the form y' = f(x,y), y(x0) = y0.
%
% y = EulerODE(f,x,y0), where
%
% f is an anonymous function representing f(x,y),
% x is a vector representing the mesh points,
% y0 is a scalar representing the initial value of y,
%
% y is the vector of solution estimates at the mesh points.
% y = 0*x;
% Pre-allocate
y(1) = y0; h = x(2)-x(1); n = length(x);
for i = 1:n-1,
    y(i+1) = y(i)+h*f(x(i),y(i));
end