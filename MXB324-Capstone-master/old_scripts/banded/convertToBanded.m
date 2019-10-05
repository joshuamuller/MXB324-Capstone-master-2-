function [ bandedA ] = convertToBanded( A )
%CONVERTTOBANDED Summary of this function goes here
%   Only works for symmetric matrices A.
    [N,~] = size(A);
	bw = getBandwidth(A);
    
    bandedA = zeros([bw,N]);

    % loop through columns
    for i = 1:N
        colSize = bw - 1;
        start = 1;
        
        if i <= bw
            colSize = i - 1;
            start = bw - colSize;
        end
        bandedA(start:start+colSize,i) = A(i-colSize:i,i);
    end
end

