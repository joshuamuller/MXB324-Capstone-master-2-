function [x, converged, k, res] = sor_csr(A, IA, JA, b, x0, tol, maxiters, w)

% Initialise
n = size(x0, 1);
k = 0;
x = x0;
normb = norm(b);
res = zeros(maxiters,1);
res(1) = norm(b - csr_mult(A, IA, JA, x)) / normb;
add_flops(flops_div);

converged = false;
[~, ii] = csr_mult(A, IA, JA, x);

while k < maxiters && res(k+1) > tol
    xold = x;

    for i = 1:n
        x(i) = b(i);
        mult1 = csr_mult(A, IA, JA, x, i, 1:i-1);
        mult2 = csr_mult(A, IA, JA, xold, i, i+1:n);
        x(i) = x(i) - mult1 - mult2;       
        x(i) = (1-w) * xold(i) + w*x(i)/ii(i);
        add_flops(2 + flops_div);
    end
    
    k = k + 1;
    res(k+1) = norm(b - csr_mult(A, IA, JA, x)) / normb;
    add_flops(flops_div);
end

if k < maxiters
    converged = true;
end

end