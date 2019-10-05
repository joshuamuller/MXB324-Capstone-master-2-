function [A] = cholesky_factorisation(A)
[n,~] = size(A);

for j = 1:n
    for i = 1:j-1
        for k = 1:i-1
            A(i,j) = A(i,j) - A(k,i) * A(k,j);
            add_flops(1); % One multiplication
        end
        A(i,j) = A(i,j) / A(i,i);
        add_flops(flops_div);
        A(j,j) = A(j,j) - A(i,j)^2;
        add_flops(flops_pow(2));
    end
    A(j,j) = sqrt(A(j,j));
    add_flops(flops_sqrt);
end