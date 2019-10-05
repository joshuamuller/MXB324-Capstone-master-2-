function [ x ] = forward_substitution_banded( A, b )
    %FORWARD_SUBSTITUTION_BANDED Summary of this function goes here
    %   Detailed explanation goes here

    [bw,N] = size(A);
    x = b;

    for j = 1:N
        if j < bw
            for i = 1:j-1
                x(j) = x(j) - A(bw-i,j)*x(j-i);
                add_flops(1);
            end
        else
            for i = 1:bw-1
                x(j) = x(j) - A(bw-i,j)*x(j-i);
                add_flops(1);
            end
        end
        
        x(j) = x(j) / A(bw,j);
        add_flops(flops_div);
    end
end

