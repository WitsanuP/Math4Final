function x = ThomasMethod(A,b)
%
% ThomasMethod uses Thomas method to find the solution vector x of a
% tridiagonal system Ax = b.
%
%   x = ThomasMethod(A,b), where
%
%       A is a tridiagonal n-by-n coefficient matrix,
%       b is the n-by-1 vector of the right-hand sides,
%
%       x is the n-by-1 solution vector.
%
n = size(A,1);
d = diag(A);        % Vector of diagonal entries of A
l = [0;diag(A,-1)]; % Vector of lower diagonal elements
u = [diag(A,1);0];  % Vector of upper diagonal elements

u(1) = u(1)/d(1); b(1) = b(1)/d(1); % First equation

for i=2:n-1,        % The next n-2 equations
    den = d(i) - u(i-1)*l(i);
    if den == 0,
        x = 'failure, division by zero';
        return
    end
 u(i) = u(i)/den; b(i) = (b(i)-b(i-1)*l(i))/den;
end

b(n)=(b(n)-b(n-1)*l(n))/(d(n)-u(n-1)*l(n)); % Last equation
x(n) = b(n);
for i=n-1:-1:1,
    x(i) = b(i) - u(i)*x(i+1);
end
x = x';