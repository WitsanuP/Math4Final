function [U, k] = PRADI(x,y,f,uleft,uright,ubottom,utop,tol,kmax)
%
% PRADI numerically solves an elliptic PDE with Dirichlet boundary
% conditions over a rectangular region using the Peaceman-Rachford
% alternating direction implicit method.
%
%   [U, k] = PRADI(x,y,f,uleft,uright,ubottom,utop,tol,kmax), where
%
%   x is the 1-by-m vector of mesh points in the x direction,
%   y is the 1-by-n vector of mesh points in the y direction,
%   f is an anonymous function representing f(x,y),
%   ubottom,uleft,utop,uright are anonymous functions describing the
%   boundary conditions,
%   tol is the tolerance (default 1e-4),
%   kmax is the maximum number of iterations (default 50),
%
%   U is the solution at the mesh points,
%   k is the number of (full) iterations needed to meet the tolerance.
%
%   Note: The default starting value at all mesh points is 0.5.
%
if nargin < 9 || isempty(kmax), kmax = 50; end
if nargin < 8 || isempty(tol), tol = 1e-4; end
y = y';
[X,Y] = meshgrid(x(2:end-1),y(2:end-1));    % Create mesh grid
m = size(X,2); n = size(X,1); N = m*n;
u = 0.5*ones(n,m);  % Starting values
h = x(2)-x(1);
% Define boundary conditions
for i = 2:m+1,
    utop_bound(i-1) = utop(x(i));
    ubottom_bound(i-1) = ubottom(x(i));
end
for i = 1:n+2,
    uleft_bound(i)=uleft(y(i));
    uright_bound(i)=uright(y(i));
end
U = [ubottom_bound;u;utop_bound]; U = [uleft_bound' U uright_bound'];
% Generate matrix A1 (first half) and A2 (second half).
A = diag(-4*ones(N,1));
v1 = diag(A,1)+1; v1(m:m:N-1) = 0;
v2 = diag(A,-1)+1; v2(n:n:N-1) = 0;
A2 = diag(v2,1)+diag(v2,-1)+A;
A1 = diag(v1,1)+diag(v1,-1)+A;
U1 = U;
for i = 1:N,
% Initialize vector b
    b0(i) = h^2*f(X(i),Y(i));
end
b0 = reshape(b0,n,m);
for k = 1:kmax,
    % First half
    b = b0-U1(1:end-2,2:end-1)-U1(3:end,2:end-1);
    b(:,1) = b(:,1)-U(2:end-1,1);
    b(:,end) = b(:,end)-U(2:end-1,end); 
    b = reshape(b',N,1);
    u = ThomasMethod(A1,b);
    % Solve tridiagonal system
    u = reshape(u,m,n);
    U1 = [U(1,2:end-1);u';U(end,2:end-1)];
    U1 = [U(:,1) U1 U(:,end)];
    % Second half
    b = b0-U1(2:end-1,1:end-2)-U1(2:end-1,3:end);
    b(1,:) = b(1,:)-U(1,2:end-1);   
    b(end,:) = b(end,:)-U(end,2:end-1);
    b = reshape(b,N,1);
    u = ThomasMethod(A2,b);
    % Solve tridiagonal system
    u = reshape(u,n,m);
    U2 = [U(1,2:end-1);u;U(end,2:end-1)];
    U2 = [U(:,1) U2 U(:,end)];
    if norm(U2-U1,inf)<=tol, break, end;
    U1 = U2;
end
    [X,Y] = meshgrid(x,y);
    U = U1;
for i = 1:n+2,
    W(i,:) = U(n-i+3,:);
    YY(i) = Y(n-i+3);
end
    U = W; Y = YY;
    surf(X,Y,U);
xlabel('x');ylabel('y');