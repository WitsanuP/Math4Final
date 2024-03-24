function u = Heat1DCN(t,x,u0,alpha)
%
% Heat1DCN numerically solves the one-dimensional heat equation, with zero
% boundary conditions, using the Crank-Nicolson method.
%
%   u = Heat1DCN(t,x,u0,alpha), where
%
%   t is the row vector of times,
%   x is the row vector of x positions,
%   u0 is the row vector of initial temperatures at the x positions,
%   alpha is a given parameter of the heat equation,
%
%   u is the approximate solution at the mesh points.
%
u = u0(:);
% u must be a column vector
k = t(2)-t(1); h = x(2)-x(1); r = (alpha/h)^2*k;

% Compute A
n = length(x);
A = diag(2*(1+r)*ones(n-2,1));
A = A + diag(diag(A,-1)-r,-1);
A = A + diag(diag(A,1)-r, 1);

% Compute B
B = diag(2*(1-r)*ones(n-2,1));
B = B + diag(diag(B,-1)+r,-1);
B = B + diag(diag(B,1)+r, 1);
C = A\B;

i = 2:length(x)-1;
for j = 1:length(t)-1,
    u(i,j+1) = C*u(i,j);
end
u = flipud(u');