function [x, converged, k] = conjugate_gradient(A, b, x0, tol, maxiters)
% CONJUGATE_GRADIENT Uses the conjugate gradient algorithm to solve the
% given system. The maxiters argument is optional; its absence will choose
% maxiters = n.

% Declarations
x = x0;
r = b - A*x;
d = r;
rTr = r' * r;
normb = norm(b);
converged = false;

if ~exist('maxiters', 'var')
    maxiters = size(A, 1);
end

% Iterate
for k = 1:maxiters
    Ad = A*d;
    alpha = rTr / (d' * Ad);
    x = x + alpha * d; 
    res = norm(b - A*x) / normb;
    if (res < tol)
        converged = true;
        break;
    end
    r = r - alpha * Ad;
    beta = 1 / rTr;
    rTr = r' * r; 
    beta = beta * rTr;
    d = r + beta * d;
end


end