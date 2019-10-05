function [x, converged, k, res] = jacobi_csr(A, IA, JA, b, x0, tol, maxiters)

% Initialise
n = size(A,1);
k = 0;
x = x0;
normb = norm(b);
res = zeros(maxiters,1);
res(1) = norm(b - csr_mult(A, IA, JA, x)) / normb;
add_flops(flops_div);

converged = false;

while k < maxiters && res(k+1) > tol
    xold = x;
    [mult, ii] = csr_mult(A, IA, JA, xold);
    x = (b - mult);
    x = x + ii.*xold;
    add_flops(size(ii, 1));
    x = x./ii;
    add_flops(flops_div * size(ii, 1));
    k = k+1;
    res(k+1) = norm(b - csr_mult(A, IA, JA, x)) / normb;
    add_flops(flops_div);
end

if k < maxiters
    converged = true;
end

end

