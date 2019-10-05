function [ bw ] = calculate_bandwidth( A, spd )
%CALCULATE_BANDWIDTH Summary of this function goes here
%   if A is square, bw is calculated.
%   if A is not square, it is assumed to be in banded (row) form.
    [n,m] = size(A);
    
    if n ~= m
        bw = n;
    else
        if nargin == 1
            spd = false;
        end

        p = 0;
        q = 0;
        for i = 1:n
            for j = n:-1:i+1
                if A(i,j) ~= 0
                    if (j-i) > q
                        q = j - i;
                    end
                end
                if A(j,i) ~= 0 && ~spd
                    if (j-i) > p
                        p = j - i;
                    end
                end
            end
        end

        bw = p + q + 1;
    end
end

