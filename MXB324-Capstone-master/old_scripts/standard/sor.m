function [x, converged, k] = sor(A, b, x0, tol, maxiters, w)

% Initialise
n = size(A,1);
k = 0;
x = x0;
normb = norm(b);
res = norm(b - A*x) / normb;
max(eig(A));

if ~exist('w', 'var')
    Tj = eye(n) - inv(diag(diag(A))) * A;
    rhoJ = max(abs(eig(Tj)));
    w = 2/(1+sqrt(1-rhoJ^2));
end

% Perform the iteration until res <= tol or maximum iterations taken
while k < maxiters && res > tol
    xold = x;
    for i = 1:n
        x(i) = b(i);
        
        x(i) = x(i) - A(i, 1:i-1)*x(1:i-1) - A(i, i+1:n)*xold(i+1:n);
        
        x(i) = (1-w)*xold(i) + w*x(i)/A(i, i);
    end
    k = k + 1;
    res = norm(b - A*x) / normb;
end

if k < maxiters
    converged = true;
else
    converged = false;
end

end