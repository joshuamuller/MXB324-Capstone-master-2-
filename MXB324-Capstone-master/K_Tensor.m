function [Kxx_n, Kxx_e, Kxx_s, Kxx_w, Kzz_n, Kzz_e, Kzz_s, Kzz_w] = K_Tensor(K_xx,K_zz, X, Z, Nx, Nz, dx, dz)
%K_Tensor calculates the Kxx and Kzz conductivities for the entire mesh
%depending on their values defined at each range of X and Z. 
% Author: Joshua

for i = 1:Nz
    for j = 1:Nx
        if X(1,j) <= 50 && Z(i,1) <= 100 && Z(i,1) >=40 %Alluvium
            Kxx(i,j) = K_xx(1);
            Kzz(i,j) = K_zz(1);
        elseif X(1,j) >= 50 && X(1,j) <= 300 && Z(i,1) <= 100 && Z(i,1) >=50
            Kxx(i,j) = K_xx(1);
            Kzz(i,j) = K_zz(1);
        end
        if X(1,j) >= 300 && X(1,j) <=500 && Z(i,1) <= 100 && Z(i,1) >=50 %Main Range Volcanics
            Kxx(i,j) = K_xx(2);
            Kzz(i,j) = K_zz(2);
        end
        if X(1,j) <= 500 && Z(i,1) <= 100 && Z(i,1) <=40 %Walloon Coal Measures
            Kxx(i,j) = K_xx(3);
            Kzz(i,j) = K_zz(3);
        elseif X(1,j) <= 500 && X(1,j)>= 300 && Z(i,1) <= 50 && Z(i,1) >=40
            Kxx(i,j) = K_xx(3);
            Kzz(i,j) = K_zz(3);
        end
        if X(1,j) >= 50 && X(1,j) <= 350 && Z(i,1) >=40 && Z(i,1) <= 50 %Confining Layer
            Kxx(i,j) = K_xx(4);
            Kzz(i,j) = K_zz(4);
        end
    end
end

Kzz_N = zeros(Nz, Nx);
Kzz_E = zeros(Nz, Nx);
Kzz_S = zeros(Nz, Nx);
Kzz_W = zeros(Nz, Nx);
Kzz_P = Kzz;
Kzz_N(2:Nz, :) = Kzz(1:Nz-1,:);
Kzz_E(:,1:Nx-1) = Kzz(:, 2:Nx);
Kzz_S(1:Nz-1,:) = Kzz(2:Nz, :);
Kzz_W(:,2:Nx) = Kzz(:,1:Nx-1);

Kxx_N = zeros(Nz, Nx);
Kxx_E = zeros(Nz, Nx);
Kxx_S = zeros(Nz, Nx);
Kxx_W = zeros(Nz, Nx);
Kxx_P = Kxx;
Kxx_N(2:Nz, :) = Kxx(1:Nz-1,:);
Kxx_E(:,1:Nx-1) = Kxx(:, 2:Nx);
Kxx_S(1:Nz-1,:) = Kxx(2:Nz, :);
Kxx_W(:,2:Nx) = Kxx(:,1:Nx-1);

Kzz_N = Kzz_N(:); Kzz_E = Kzz_E(:); Kzz_S = Kzz_S(:); Kzz_W = Kzz_W(:); Kzz_P = Kzz_P(:);
Kxx_N = Kxx_N(:); Kxx_E = Kxx_E(:); Kxx_S = Kxx_S(:); Kxx_W = Kxx_W(:); Kxx_P = Kxx_P(:);

for i = 1:length(Kxx_N)
    Kzz_n(i) = dx(2,2)/2 *(Kzz_N(i) + Kzz_P(i));
    Kzz_e(i) = dx(2,2)/2 *(Kzz_E(i) + Kzz_P(i));
    Kzz_s(i) = dx(2,2)/2 *(Kzz_S(i) + Kzz_P(i));
    Kzz_w(i) = dx(2,2)/2 *(Kzz_W(i) + Kzz_P(i));
    Kxx_n(i) = dz(2,2)/2 *(Kxx_N(i) + Kxx_P(i));
    Kxx_e(i) = dz(2,2)/2 *(Kxx_E(i) + Kxx_P(i));
    Kxx_s(i) = dz(2,2)/2 *(Kxx_S(i) + Kxx_P(i));
    Kxx_w(i) = dz(2,2)/2 *(Kxx_W(i) + Kxx_P(i));
end
end