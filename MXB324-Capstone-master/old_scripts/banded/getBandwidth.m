function [ bw ] = getBandwidth( A )
[N,~] = size(A);
%GETBANDWIDTH Calculates the bandwidth of a square matrix matrix A.
    bw = 0;
    for i = 1:N
        for j = N:-1:i
            if A(i,j) ~= 0
                rowSize = j-i + 1;
                if rowSize > bw
                    bw = rowSize;
                end
            end
        end
    end
end

