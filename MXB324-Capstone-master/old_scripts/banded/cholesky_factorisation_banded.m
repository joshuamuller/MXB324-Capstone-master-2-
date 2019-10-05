function [ A ] = cholesky_factorisation_banded( A )
    [bw,N] = size(A);
    
    for j = 1:N
        % in first cases, have to only start with values in band
        if j < bw
            first = bw - j;
            for i = 1:j-1
                mapped_i = first + i;
                
                for k = 1:i-1
                    mapped_k = first + k;
                    offset = j - i;
                    A(mapped_i,j) = A(mapped_i,j) - A(mapped_k + offset,i) * A(mapped_k,j);
                    add_flops(1);
                end
                A(mapped_i, j) = A(mapped_i, j) / A(bw, i);
                add_flops(flops_div);
                A(bw,j) = A(bw,j) - A(mapped_i,j)^2;
                add_flops(flops_pow(2));
            end
        else
            for i = 1:bw-1
                mapped_i = j - bw + i;
                for k = 1:i-1
                    offset = bw - i;
                    A(i,j) = A(i,j) - A(k + offset,mapped_i) * A(k,j);
                    add_flops(1);
                end
               
                A(i, j) = A(i, j) / A(bw, mapped_i);
                add_flops(flops_div);
                A(bw,j) = A(bw,j) - A(i,j)^2;
                add_flops(flops_pow(2));
            end
        end
        A(bw,j) = sqrt(A(bw,j));
        add_flops(flops_sqrt);
    end

end

