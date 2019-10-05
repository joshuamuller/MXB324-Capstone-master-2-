function X = linear_solve_spd(A, b)

Aout = cholesky_factorisation(A);

G = triu(Aout)';

X2 = forward_substitution(G,b);
X = backward_substitution(G', X2);