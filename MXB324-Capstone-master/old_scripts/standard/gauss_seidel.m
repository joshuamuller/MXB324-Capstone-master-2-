function [x, converged, k] = gauss_seidel(A, b, x0, tol, maxiters)
% GAUSS_SEIDEL Gauss-Seidel method
% [x, converged, k] = gauss_seidel(A, b, x0, tol, maxiters) performs
% Gauss-Seidel method iteration to solve A*x=b, starting with x = x0 and 
% iterating until either the number of iterations k equals maxiters, or 
% until norm(b-A*x)/norm(b) <= tol. The value of converged is true if the 
% method converged, or false otherwise.

% Initialise
n = size(A,1);
k = 0;
x = x0;
normb = norm(b);
converged = false;
res = norm(b - A*x) / normb;

% Perform the iteration until res <= tol or maximum iterations taken
while k < maxiters && res > tol
    xold = x;
    for i = 1:n
        x(i) = b(i);
        
        x(i) = x(i) - A(i, 1:i-1)*x(1:i-1) - A(i, i+1:n)*xold(i+1:n);
        
        x(i) = x(i)/A(i, i);
    end
    k = k + 1;
    res = norm(b - csr_mult(A, IA, JA, x)) / normb;
end

if k < maxiters
    converged = true;
end

end