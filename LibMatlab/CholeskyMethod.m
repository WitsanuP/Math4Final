function x = CholeskyMethod(A,b)
%
% CholeskyMethod uses the Cholesky factorization of matrix A and
% solves the ensuing triangular systems to find the solution vector x.
%
%   x = CholeskyMethod(A,b), where
%
%       A is a symmetric, positive definite n-by-n coefficient matrix,
%       b is the n-by-1 vector of the right-hand sides,
%
%       x is the n-by-1 solution vector.
%
[L, U] = CholeskyFactor(A);     % Find Cholesky factorization of A
n = size(A,1);

% Solve the lower triangular system Ly = b (forward substitution)
y = zeros(n,1);
y(1) = b(1)/L(1,1);
for i = 2:n,
    y(i) = (b(i)-L(i,1:i-1)*y(1:i-1))/L(i,i);
end

% Solve the upper triangular system L'x = y (back substitution)
x = zeros(n,1);
x(n) = y(n)/U(n,n);
    for i = n-1:-1:1,
        x(i) = (y(i)-U(i,i+1:n)*x(i+1:n))/U(i,i);
    end
end