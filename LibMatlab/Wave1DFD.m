function u = Wave1DFD(t,x,u0,ut0,alpha)
%
% Wave1DFD numerically solves the one-dimensional wave equation, with zero
% boundary conditions, using the finite-difference method.
%
%   u = Wave1DFD(t,x,u0,ut0,alpha), where
%
%   t is the row vector of times,
%   x is the row vector of x positions,
%   u0 is the row vector of initial displacements for x positions,
%   ut0 is the row vector of initial velocities for x positions,
%   alpha is a given parameter of the wave equation,
%
%   u is the approximate solution at the mesh points.
%
u = u0(:); ut = ut0(:); % u and ut must be column vectors
k = t(2)-t(1); h = x(2)-x(1); r = (k*alpha/h)^2;
if r > 1,
    warning('Method is unstable and divergent. Results will be inaccurate.')
end
i = 2:length(x)-1;
u(i,2) = (1-r)*u(i,1) + r/2*(u(i-1,1) + u(i+1,1)) + k*ut(i);
for j = 2:length(t)-1,
    u(i,j+1) = -u(i,j-1) + 2*(1-r)*u(i,j) + r*(u(i-1,j) + u(i+1,j));
end
u = flipud(u');