function X = linear_solve_spd_packed(A, b)

A = cholesky_factorisation_packed(A);

X = forward_substitution_packed(A,b);
X = backward_substitution_packed(A, X);