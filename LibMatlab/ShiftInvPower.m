function [e_val, e_vec, k] = ShiftInvPower(A, alpha, x1, tol, kmax)
%
% ShiftInvPower uses the shifted inverse power method to find the
% eigenvalue of a matrix closest to a specified value. It also returns
% the eigenvector associated with this eigenvalue.
%
%   [e_val, e_vec, k] = ShiftInvPower(A, alpha, x1, tol, kmax), where
%
%   A is an n-by-n matrix,
%   alpha is a specied value,
%   x1 is the n-by-1 initial vector (default ones),
%   tol is the tolerance (default 1e-4),
%   kmax is the maximum number of iterations (default 50),
%
%   e_val is the estimated eigenvalue,
%   e_vec is the corresponding eigenvector,
%   k is the number of iterations required for convergence.
%
n = size(A,1);
if nargin < 3 || isempty(x1), x1 = ones(n,1); end
if nargin < 4 || isempty(tol), tol = 1e-4; end
if nargin < 5 || isempty(kmax), kmax = 50; end
x(:,1) = x1./norm(x1,2);
betas(1) = 0;
for k = 1:kmax,
    x(:,k+1) = DoolittleMethod(A-alpha*eye(n,n),x(:,k));    
    betas(k+1) = x(:,k).'*x(:,k+1);
    x(:,k+1) = x(:,k+1)./norm(x(:,k+1),2);
    if abs(betas(k+1)-betas(k)) < tol,  % Check for convergence
        break
    end
end
betas = betas(end);
e_val = 1/betas+alpha;
e_vec = x(:, end);