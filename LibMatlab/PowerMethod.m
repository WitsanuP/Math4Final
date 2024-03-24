function [e_val, e_vec, k] = PowerMethod(A, x1, tol, kmax)
%
% PowerMethod approximates the dominant eigenvalue and the corresponding
% eigenvector of a square matrix.
%
% [e_val, e_vec, k] = PowerMethod(A, x1, tol, kmax), where
%
%   A is an n-by-n matrix,
%   x1 is the n-by-1 initial vector (default ones),
%   tol is the tolerance (default 1e-4),
%   kmax is the maximum number of iterations (default 50),
%
%   e_val is the approximated dominant eigenvalue,
%   e_vec is the corresponding eigenvector,
%   k is the number of iterations required for convergence.
%
n = size(A,1);
if nargin < 2 || isempty(x1), x1 = ones(n,1); end
if nargin < 3 || isempty(tol), tol = 1e-4; end
if nargin < 4 || isempty(kmax), kmax = 50; end
x(:,1) = x1./norm(x1, 2);
x(:,2) = A*x(:,1); alpha(2) = x(:,1).'*x(:,2);
x(:,2) = x(:,2)./norm(x(:,2),2);

for k = 2:kmax,
    x(:,k+1) = A*x(:,k); % Generate next vector
    alpha(k+1) = x(:,k).'*x(:,k+1);
    x(:,k+1) = x(:,k+1)./norm(x(:,k+1),2);
    if abs(alpha(k+1)-alpha(k)) < tol, % Terminating condition
        break
    end
end
e_val = alpha(end); e_vec = x(:,end);