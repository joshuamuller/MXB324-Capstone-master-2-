function [ X ] = linear_solve_spd_banded( A, b)
    A = cholesky_factorisation_banded(A);
    X = forward_substitution_banded(A,b);
    X = backward_substitution_banded(A, X);
end

