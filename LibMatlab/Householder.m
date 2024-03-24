function T = Householder(A)
%
% Householder uses Householder's method to transform a symmetric matrix
% into a tridiagonal matrix.
%
%   T = Householder(A), where
%
%   A is an n-by-n real, symmetric matrix,
%
%   T is an n-by-n tridiagonal matrix.
%
N = size(A,1);
for n = 1:N-2,
    S = sqrt(A(n+1:end,n)'*A(n+1:end,n)); % Compute sigma
    v(1:N,1) = 0; % Set initial set of entries to 0
    v(n+1) = sqrt((1+abs(A(n+1,n))/S)/2); % First non-zero entry
    sn = sign(A(n+1,n)); % Determine sign of relevant entry
    v(n+2:N) = sn*A(n+2:end,n)/2/v(n+1)/S; % Compute remaining entries
    P = eye(N)-2*(v*v'); % Compute the symmetric, orthogonal matrices
    A = P\A*P; % Compute sequence of matrices
end
T = A;