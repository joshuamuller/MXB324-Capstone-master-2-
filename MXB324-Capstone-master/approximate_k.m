function [k_n, k_s, k_e, k_w] = approximate_k(H, k, avg_type, dz, dx, Nz, Nx)

k = reshape(k,[Nz,Nx]);
k_N = zeros(Nz, Nx);
k_E = zeros(Nz, Nx);
k_S = zeros(Nz, Nx);
k_W = zeros(Nz, Nx);
k_P = k;
k_N(2:Nz, :) = k(1:Nz-1,:);
k_E(:,1:Nx-1) = k(:, 2:Nx);
k_S(1:Nz-1,:) = k(2:Nz, :);
k_W(:,2:Nx) = k(:,1:Nx-1);

k_N = k_N(:); k_E = k_E(:); k_S = k_S(:); k_W = k_W(:); k_P = k_P(:);

%
H = reshape(H,[Nz,Nx]);
H_N = zeros(Nz, Nx);
H_E = zeros(Nz, Nx);
H_S = zeros(Nz, Nx);
H_W = zeros(Nz, Nx);
H_P = H;
H_N(2:Nz, :) = H(1:Nz-1,:);
H_E(:,1:Nx-1) = H(:, 2:Nx);
H_S(1:Nz-1,:) = H(2:Nz, :);
H_W(:,2:Nx) = H(:,1:Nx-1);

H_N = H_N(:); H_E = H_E(:); H_S = H_S(:); H_W = H_W(:); H_P = H_P(:); H = H(:);

%%
for i = 1:length(H)
if avg_type == "arithmetic"
    
    k_n(i) = (k_P(i) + k_N(i))/2;
    k_s(i) = (k_P(i) + k_S(i))/2;
    k_e(i) = (k_P(i) + k_E(i))/2;
    k_w(i) = (k_P(i) + k_W(i))/2;

elseif avg_type == "harmonic"
    
    k_n = dz(2,2)*(dz(2,2)/(2*k_N(i)) + dz(2,2)/(2*k_P(i)))^-1;
    k_s = dz(2,2)*(dz(2,2)/(2*k_S(i)) + dz(2,2)/(2*k_P(i)))^-1;
    k_e = dx(2,2)*(dx(2,2)/(2*k_E(i)) + dx(2,2)/(2*k_P(i)))^-1;
    k_w = dx(2,2)*(dx(2,2)/(2*k_W(i)) + dx(2,2)/(2*k_P(i)))^-1;
    
elseif avg_type == "upstream"
    
    if H_P(i) < H_W(i)
        k_w(i) = k_W(i);
    else
        k_w(i) = k_P(i);
    end
    
    if H_P(i) < H_E(i)
        k_e(i) = k_E(i);
    else
        k_e(i) = k_P(i);
    end
        
    if H_P(i) < H_N(i)
        k_n(i) = k_N(i);
    else
        k_n(i) = k_P(i);
    end
    
    if H_P(i) < H_S(i)
        k_s(i) = k_S(i);
    else
        k_s(i) = k_P(i);
    end 
    
end
end
end