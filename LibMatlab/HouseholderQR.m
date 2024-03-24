function [T, Tfinal, e_vals, m] = HouseholderQR(A, tol, kmax)
%
% HouseholderQR uses Householder's method and repeated applications of
% QR factorization to estimate the eigenvalues of a real, symmetric matrix.
%
%   [T, Tfinal, e_vals, m] = HouseholderQR(A, tol, kmax), where
%
%   A is an n-by-n real, symmetric matrix,
%   tol is the tolerance used in the QR process (default 1e-4),
%   kmax is the maximum number of QR iterations (default 50),
%
%   T is the tridiagonal matrix created by Householder's method,
%   Tfinal is the final tridiagonal matrix,
%   e_vals is a list of estimated eigenvalues of matrix A,
%   m is the number of iterations needed for convergence.
%
% Note that this function calls the user-defined function Householder!
%
if nargin < 2 || isempty(tol), tol = 1e-4; end
if nargin < 3 || isempty(kmax), kmax = 50; end
T = Householder(A); % Call Householder
T(:,:,1) = T;
% QR factorization to reduce the off-diagonal entries of T
for m = 1:kmax,
    [Q,R] = qr(T(:,:,m));
    T(:,:,m+1) = R*Q;
    % Compare diagonals of two successive T matrices
    if norm(diag(T(:,:,m+1))-diag(T(:,:,m))) < tol,
        break;
    end
end
Tfinal = T(:,:,end); T = T(:,:,1); e_vals = diag(Tfinal);