function A = cholesky_factorisation_packed(A)
    [n,~] = size(A);
    n = -(1 - sqrt(1 + n*8))*0.5;
    add_flops(2 + flops_sqrt);
    for j = 1:n
        mappedAjj = apply_mapping(j, j);
        for i = 1:j-1
            mappedAij = apply_mapping(i, j);
            mappedAii = apply_mapping(i, i);
            for k = 1:i-1
                A(mappedAij) = A(mappedAij) - A(apply_mapping(k,i)) * A(apply_mapping(k,j));
                add_flops(1);
            end
            
            A(mappedAij) = A(mappedAij) / A(mappedAii);
            A(mappedAjj) = A(mappedAjj) - A(mappedAij)^2;
            add_flops(flops_div + flops_pow(2));
        end
        A(mappedAjj) = sqrt(A(mappedAjj));
        add_flops(flops_sqrt);
    end
end 

function idx = apply_mapping(i,j)
    idx = i + j*(j-1)*0.5;
    add_flops(2); % Multiplication and division.
end