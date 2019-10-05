Nx = 51; Nz = 21;

[SM, X, Z, Mesh_ref, H_P, H_N, H_S, H_W, H_E] = Mesh_creator(Nx,Nz);

K_xx = [2.6 0.08 3.9 4.3]; %meters per day
K_zz = [0.91 0.016 1.17 0.21]; %meters per day
Psi_res = [0.01 0.106 0.01 0.01];
Psi_sat = [0.33 0.4686 0.1 0.075];
alpha   = [1.43 1.04 2.8 2.5];
n       = [1.51 1.3954 2.239 2.0];

K = [K_xx K_zz];