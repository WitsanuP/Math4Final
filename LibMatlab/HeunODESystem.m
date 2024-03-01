function u = HeunODESystem(f,x,u0)
%
% HeunODESystem uses Heun's method to solve a system of first-order
% initial-value problems in the form u' = f(x,u), u(x0) = u0.
%
% u = HeunODESystem(f,x,u0), where
%
%   f is an anonymous m-dim. vector function representing f(x,u),
%   x is an (n+1)-dim. vector representing the mesh points,
%   u0 is an m-dim. vector representing the initial state vector,
%
%   u is an m-by-(n+1) matrix, each column the vector of solution
%   estimates at a mesh point.
%
u(:,1) = u0;    % The first column is set to be the initial vector u0
h = x(2) - x(1); n = length(x);
for i = 1:n-1,
     k1 = f(x(i),u(:,i));
     k2 = f(x(i)+h,u(:,i)+h*k1);
    u(:,i+1) = u(:,i)+h*(k1+k2)/2;
end