function [X, Z, Mesh_ref, N_P, N_N, N_S, N_W, N_E] = Mesh_creator(Nx,Nz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Nodes = Nx*Nz;
Numbering = 1:Nx*Nz;
Mesh_ref = reshape(Numbering, [Nz, Nx]); %Reference Mesh
[X,Z] = meshgrid(linspace(0,500,Nx),linspace(0,100,Nz));
X = flipud(X); Z = flipud(Z);
N_P = Mesh_ref(2:Nz-1, 2:Nx-1);
N_P = N_P(:);
N_N = Mesh_ref(1:Nz-2, 2:Nx-1);
N_N = N_N(:);
N_S = Mesh_ref(3:Nz, 2:Nx-1);
N_S = N_S(:);
N_W  = Mesh_ref(2:Nz-1, 1:Nx-2);
N_W = N_W(:);
N_E  = Mesh_ref(2:Nz-1, 3:Nx);
N_E = N_E(:);
end

