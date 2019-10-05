% Evaluate the flux (q) at north, south, east and west control volume
% faces. This function returns four vectors q_n, q_s, q_e, q_w which are
% the same length as the number of nodes in the mesh. q_e(i) represents the
% flux through the east face for node i.  On boundaries, the q's represent 
% the boundary flux i.e. q = 0 for no flow, q = q_rain for rainfall.

function [q_n, q_s, q_e, q_w] = evaluate_q(head_profile)

%% Initialise arrays to store fluxes

q_n = zeros(1, Nx*Nz);
q_s = zeros(1, Nx*Nz);
q_e = zeros(1, Nx*Nz);
q_w = zeros(1, Nx*Nz);

%% Internal Nodes

for i = 1:length(N_P)

    q_n(i) = -1/2 * k_n(i) * Kzz_n(i) * (H(N_N(i)) - H(N_P(i)))/dz;
    q_s(i) = -1/2 * k_s(i) * Kzz_s(i) * (H(N_P(i)) - H(N_S(i)))/dz;
    q_e(i) = -1/2 * k_e(i) * Kxx_e(i) * (H(N_E(i)) - H(N_P(i)))/dx;
    q_w(i) = -1/2 * k_w(i) * Kxx_w(i) * (H(N_P(i)) - H(N_W(i)))/dx;

end

%% Bottom Left Corner (0, 0)

q_n(Nz) = -1/2 * k_n(Nz) * Kzz_n(Nz) * (H(Nz-1) - H(Nz))/dz;
q_s(Nz) = 0;
q_e(Nz) = -1/2 * k_e(Nz) * Kxx_e(Nz) * (H(2*Nz) - H(Nz))/dx;
q_w(Nz) = 0;

%% Bottom Right Corner (500, 0)

n = Nx*Nz;

q_n(n) = -1/2 * k_n(n) * Kzz_n(n) * (H(n-1) - H(n))/dz;
q_s(n) = 0;
q_w(n) = -1/2 * k_w(n) * Kxx_w(n) * (H(n) - H(n-Nz))/dx;

if CSG == true
    dH_CSG = (20 - H(n))/5000;
    q_e(n) = 3.9 * dH_CSG; % CSG extraction
else
    q_e(n) = 0; % rock
end

%% Top Left Corner (0, 100)

q_s(1) = -1/2 * k_s(1) * Kzz_s(1) * (H(1) - H(2))/dz;
q_e(1) = -1/2 * k_e(1) * Kxx_e(1) * (H(1+Nz) - H(1))/dx;

if S(1) < 1
    q_n(1) = q_rain; % rain enters soil
elseif S(1) == 1
    q_n(1) = 0; % soil already saturated
else
    disp("Soil saturation exceeds 1");
    q_n(1) = NaN;
end

if river == true
    % put river condition here
else
    q_w(1) = 0; % no river
end

%% Top Right Corner (500, 100)

q_s(Mesh_ref(1, end)) = -1/2 * k_s(Mesh_ref(1, end)) * Kzz_s(Mesh_ref(1, end)) * (H(Mesh_ref(1, end)) - H(Mesh_ref(2, end)))/dz;
q_w(Mesh_ref(1, end)) = -1/2 * k_w(Mesh_ref(1, end)) * Kxx_w(Mesh_ref(1, end)) * (H(Mesh_ref(1, end)) - H(Mesh_ref(1, end-1)))/dx;

q_e(Mesh_ref(1, end)) = 0; % rock

if S(Mesh_ref(1, end)) < 1
    q_n(Mesh_ref(1, end)) = q_rain; % rain enters soil
elseif S(Mesh_ref(1, end)) == 1
    q_n(Mesh_ref(1, end)) = 0; % soil already saturated
else
    disp("Soil saturation exceeds 1");
    q_n(Mesh_ref(1, end)) = NaN;
end

%% Left boundary (x = 0)

for i = 2:Nz-1
    
    q_n(i) = -1/2 * k_n(i) * Kzz_n(i) * (H(i-1) - H(i))/dz;
    q_s(i) = -1/2 * k_s(i) * Kzz_s(i) * (H(i) - H(i+1))/dz;
    q_e(i) = -1/2 * k_e(i) * Kxx_e(i) * (H(i+Nz) - H(i))/dx;
    
    if Z(i) > 75 && river == true
        % put river condition here;
    else
        q_w(i) = 0; % rock
    end
          
end

%% Right boundary (x = 500)

for i = 2:Nz-1
    
    q_n(Mesh_ref(i, end)) = -1/2 * k_n(Mesh_ref(i, end)) * Kzz_n(Mesh_ref(i, end)) * (H(Mesh_ref(i-1, end)) - H(Mesh_ref(i, end)))/dz;
    q_s(Mesh_ref(i, end)) = -1/2 * k_s(Mesh_ref(i, end)) * Kzz_s(Mesh_ref(i, end)) * (H(Mesh_ref(i, end)) - H(Mesh_ref(i+1, end)))/dz;
    q_w(Mesh_ref(i, end)) = -1/2 * k_w(Mesh_ref(i, end)) * Kxx_w(Mesh_ref(i, end)) * (H(Mesh_ref(i, end)) - H(Mesh_ref(i, end-1)))/dx;
    
    if Z(i) < 5 && CSG == true
        dH_CSG = (20 - H(Mesh_ref(i, end)))/5000;
        q_e(Mesh_ref(i, end)) = 3.9 * dH_CSG; % CSG extraction
    else
        q_e(Mesh_ref(i, end)) = 0; % rock
    end
          
end

%% Bottom boundary (z = 0)

for i = 2:Nx-1
    
    q_n(Mesh_ref(end, i)) = -1/2 * k_n(Mesh_ref(end, i)) * Kzz_n(Mesh_ref(end, i)) * (H(Mesh_ref(end-1, i)) - H(Mesh_ref(end, i)))/dz;
    q_s(Mesh_ref(end, i)) = 0; % rock
    q_e(Mesh_ref(end, i)) = -1/2 * k_e(Mesh_ref(end, i)) * Kxx_e(Mesh_ref(end, i)) * (H(Mesh_ref(end, i+1)) - H(Mesh_ref(end, i)))/dx;
    q_w(Mesh_ref(end, i)) = -1/2 * k_w(Mesh_ref(end, i)) * Kxx_w(Mesh_ref(end, i)) * (H(Mesh_ref(end, i)) - H(Mesh_ref(end, i-1)))/dx;

end
%% Top boundary (z = 100)
for i = 2:Nx-1

    q_s(Mesh_ref(1, i)) = -1/2 * k_s(Mesh_ref(1, i)) * Kzz_s(Mesh_ref(1, i)) * (H(Mesh_ref(1, i)) - H(Mesh_ref(2, i)))/dz;
    q_e(Mesh_ref(1, i)) = -1/2 * k_e(Mesh_ref(1, i)) * Kxx_e(Mesh_ref(1, i)) * (H(Mesh_ref(1, i-1)) - H(Mesh_ref(1, i)))/dx;
    q_w(Mesh_ref(1, i)) = -1/2 * k_w(Mesh_ref(1, i)) * Kxx_w(Mesh_ref(1, i)) * (H(Mesh_ref(1, i)) - H(Mesh_ref(1, i+1)))/dx;

    if S(Mesh_ref(1, i)) < 1
        q_n(Mesh_ref(1, i)) = q_rain; % rain enters soil
    elseif S(Mesh_ref(1, i)) == 1
        q_n(Mesh_ref(1, i)) = 0; % soil already saturated
    else
        disp("Soil saturation exceeds 1");
        q_n(Mesh_ref(1, i)) = NaN;
    end
end
end
