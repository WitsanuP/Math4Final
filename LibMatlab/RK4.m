function y = RK4(f,x,y0)
%
% RK4 uses the classical RK4 method to solve a first-order initial-value
% problem in the form y' = f(x,y), y(x0) = y0.
%
% y = RK4(f,x,y0), where
%
% f is an anonymous function representing f(x,y),
% x is a vector representing the mesh points,
% y0 is a scalar representing the initial value of y,
%
% y is the vector of solution estimates at the mesh points.
%
y = 0*x; % Pre-allocate
y(1) = y0; h = x(2)-x(1); n = length(x);
for i = 1:n-1,
    k1 = f(x(i),y(i));
    k2 = f(x(i)+h/2,y(i)+h*k1/2);
    k3 = f(x(i)+h/2,y(i)+h*k2/2);
    k4 = f(x(i)+h,y(i)+h*k3);
    y(i+1) = y(i)+h*(k1+2*k2+2*k3+k4)/6;
end