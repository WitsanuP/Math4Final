function x = GaussPivotScale(A,b)
%
% GaussPivotScale uses Gauss elimination with partial pivoting and
% row scaling to solve the linear system Ax = b.
%
%   x = GaussPivotScale(A,b), where
%
%   A is the n-by-n coefficient matrix,
%   b is the n-by-1 result vector,
%
%   x is the n-by-1 solution vector.
%
n = length(b);
A = [A b];      % Define augmented matrix
for k = 1:n-1,
    % Find maximum magnitude in each row
    M = max(abs(A(k:end, k:end-1)), [], 2);
    a = abs(A(k:end, k));   % Find maximum in kth column
    I = max(a./M);  % Find row with maximum ratio
    I = I + k - 1;  % Adjust relative row to actual row
    if I > k
        A([k I], :)= A([I k], :);   % Pivot rows
    end
    m = A(k+1:n, k)/A(k, k);    % Construct multipliers
    [Ak, M] = meshgrid(A(k, :), m); % Create mesh
    A(k+1:n, :) = A(k+1:n, :) - Ak.*M;
end
Ab = A;
% Find the solution vector using back substitution
x = BackSub(Ab);