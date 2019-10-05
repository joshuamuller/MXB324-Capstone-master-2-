function [x, converged, k, res] = conjugate_gradient_csr(A, IA, JA, b, x0, tol, maxiters)
% CONJUGATE_GRADIENT Uses the conjugate gradient algorithm to solve the
% given system. The maxiters argument is optional; its absence will choose
% maxiters = n.

% Declarations
x = x0;
r = b - csr_mult(A, IA, JA, x);
d = r;
rTr = r' * r;
add_flops(size(r, 1));
normb = norm(b);
converged = false;

if ~exist('maxiters', 'var')
    maxiters = size(A, 1);
end

res = zeros(maxiters,1);

% Iterate
for k = 1:maxiters
    Ad = csr_mult(A, IA, JA, d);
    alpha = rTr / (d' * Ad);
    add_flops(flops_div + size(Ad, 1));
    x = x + alpha * d; 
    add_flops(1);
    res(k) = norm(b - csr_mult(A, IA, JA, x)) / normb;
    add_flops(flops_div);
    if (res(k) < tol)
        converged = true;
        break;
    end
    r = r - alpha * Ad;
    add_flops(size(Ad, 1));
    beta = 1 / rTr;
    add_flops(flops_div);
    rTr = r' * r; 
    add_flops(size(r, 1));
    beta = beta * rTr;
    d = r + beta * d;
    add_flops(2);
end


end