function [x, k, MGSnorm] = GaussSeidel(A, b, x0, tol, kmax)
% GaussSeidel uses the Gauss-Seidel iteration method to approximate the solution of Ax = b.
%   [x, k, MGSnorm] = GaussSeidel (A, b, x0, tol, kmax), where
%
%   A is the n-by-n coefficient matrix,
%   b is the n-by-1 right-hand side vector,
%   x0 is the n-by-1 initial vector (default zeros),
%   tol is the scalar tolerance for convergence (default 1e-4),
%
%   kmax is the maximum number of iterations (default 100), 
%   x is the n-by-1 solution vector,
%   k is the number of iterations required for convergence,
%   MGSnorm is the infinite norm of the Gauss-Seidel iteration matrix.
%
if nargin < 3 || isempty(x0), x0 = zeros(size(b)); end
if nargin < 4 || isempty(tol), tol = 1e-4; end
if nargin < 5 || isempty(kmax), kmax = 100; end
x(:, 1) = x0;
D = diag(diag(A)); At = A - D;
L = tril(At); U = At - L;
% Norm of Gauss-Seidel iteration matrix
M = -(D + L)\U; MGSnorm = norm(M, inf); B = (D + L)\b;
% Perform iterations up to kmax

for k = 1:kmax,
    x(:, k+1) = M*x(:, k) + B;
    if norm(x(:,k+1)-x(:, k)) < tol,
        break
    end
end
x = x(:, end);  