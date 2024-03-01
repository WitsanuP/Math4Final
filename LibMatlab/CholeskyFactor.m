function [L, U] = CholeskyFactor(A)
%
% CholeskyFactor returns the Cholesky factorization of matrix A.
%
%   [L, U] = CholeskyFactor(A), where
%
%       A is a symmetric, positive definite n-by-n matrix,
%
%       L is a lower triangular matrix,
%       U = L' is an upper triangular matrix.
%
n = size(A,1);
L = zeros(n,n); % Initialize
for i = 1:n,
    L(i,i) = sqrt(A(i,i)-L(i,1:i-1)*L(i,1:i-1)');
    for j = i+1:n,
        L(j,i) = (A(j,i)-L(j,1:i-1)*L(i,1:i-1)')/L(i,i);
    end
end
U = L';