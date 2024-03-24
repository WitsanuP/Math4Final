function u = Heat1DFD(t,x,u0,alpha)
%
% Heat1DFD numerically solves the one-dimensional heat equation, with zero
% boundary conditions, using the finite-difference method.
%
%   u = Heat1DFD(t,x,u0,alpha), where
%
%   t is the row vector of times,
%   x is the row vector of x positions,
%   u0 is the row vector of initial temperatures at the x positions,
%   alpha is a given parameter of the heat equation,
%
%   u is the approximate solution at the mesh points.
%
u = u0(:);  % u must be a column vector
k = t(2)-t(1); h = x(2)-x(1); r = (alpha/h)^2*k;
if r > 0.5,
    warning('Method is unstable and divergent. Results will be inaccurate.')
end
i = 2:length(x)-1;
for j = 1:length(t)-1,
    u(i,j+1) = (1-2*r)*u(i,j) + r*(u(i-1,j)+u(i+1,j));
end
u = flipud(u');