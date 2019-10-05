function [ x ] = backward_substitution_banded( A, b )
%BACKWARD_SUBSTITUTION_BANDED Summary of this function goes here
%   Detailed explanation goes here

	[bw,N] = size(A);
    x = b;

    for j = N:-1:1
        % Divide current
        x(j) = x(j) / A(bw,j);
        add_flops(flops_div);
        
        % Divide previous rows
        last = bw-1;
        if j < bw
            last = j-1;
        end
        for i = 1:last
        	x(j-i) = x(j-i) - A(bw-i,j)*x(j);
            add_flops(1);
        end
    end
end

