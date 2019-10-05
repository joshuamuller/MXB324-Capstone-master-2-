function [ p_out ] = reverseRCM( p )
%REVERSERCM Summary of this function goes here
%   Detailed explanation goes here
    n = size(p,2);
    p_out = zeros(size(p));
    
    for i = 1:n
        p_out(i) = find(p==i);
    end
end

