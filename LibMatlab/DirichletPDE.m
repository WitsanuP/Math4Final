function U = DirichletPDE(x,y,f,uleft,uright,ubottom,utop)
%
%
%   DirichletPDE numerically solves an elliptic PDE with Dirichlet boundary
%   conditions over a rectangular region.
%
%   U = DirichletPDE(x,y,f,uleft,uright,ubottom,utop), where
%
%   x is the 1-by-m vector of mesh points in the x direction,
%   y is the 1-by-n vector of mesh points in the y direction,
%   f is an anonymous function representing f(x,y),
%   ubottom(x),utop(x),uright(y),uleft(y) are anonymous functions
%   describing the boundary conditions,
%
%   U is the solution at the interior mesh points.
%
y = y'; m = size(x,2);n = size(y,1); N = (m-2)*(n-2);
A = diag(-4*ones(N,1));
A = A + diag(diag(A,n-2)+1,n-2);
A = A + diag(diag(A,2-n)+1,2-n);
v = ones(N-1,1);    % Create vector of ones
v(n-2:n-2:end) = 0; % Insert zeros
A = A + diag(v,1);  % Add upper diagonal
A = A + diag(v,-1); % Add lower diagonal
[X,Y] = meshgrid(x(2:end-1),y(end-1:-1:2)); % Create mesh
h = x(2)-x(1);
% Define boundary conditions
for i = 2:m-1,
    utop_bound(i-1) = utop(x(i));
    ubottom_bound(i-1) = ubottom(x(i));
end
for i = 1:n,
    uleft_bound(i) = uleft(y(n+1-i));
    uright_bound(i) = uright(y(n+1-i));
end
b = 0;
% Initialize vector b
for i = 1:N,
    b(i) = h^2*f(X(i),Y(i));
end
b(1:n-2:N) = b(1:n-2:N)-utop_bound;
b(n-2:n-2:N) = b(n-2:n-2:N)-ubottom_bound;
b(1:n-2) = b(1:n-2)-uleft_bound(2:n-1);
b(N-(n-3):N) = b(N-n+3:N)-uright_bound(2:n-1);
u = A\b';   % Solve the system
U = reshape(u,n-2,m-2);
U = [utop_bound;U;ubottom_bound];
U = [uleft_bound' U uright_bound'];
[X,Y] = meshgrid(x,y(end:-1:1));
surf(X,Y,U);    % 3D plot of the numerical results
xlabel('x');ylabel('y');