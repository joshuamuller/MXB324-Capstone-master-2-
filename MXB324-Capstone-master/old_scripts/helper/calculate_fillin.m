function [ fill_in ] = calculate_fillin( A, A_chol )
%CALCULATE_FILLIN Summary of this function goes here
%   Detailed explanation goes here
    mask = A_chol ~= 0 & ~( (A_chol ~= 0) & (A ~= 0));
    
    fill_in = nnz(mask);
end

