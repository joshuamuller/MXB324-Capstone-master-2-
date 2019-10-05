function [F, q_n, q_e, q_s, q_w] = system_function(h, N, Mesh_ref, N_P, ...
    N_N, N_S, N_W, N_E, dx, dz, L, Ks, bore, river, vegetation, ...
    rain_percentages, rainfall, S, psi_old, psi, q_n_old, q_e_old, ...
    q_s_old, q_w_old, ns, ms, alpha, t)
%SYSTEM_FUNCTION Evaluates the deviation from conservation for the given h.
%
%Task name: system_function
%Authors: Joshua, Matthew, Paul
%
%Revisions [Date, people, reason]:
%- 21/09/2019, Paul, Reworked for proper system function and Jacobian FD.

F = zeros(length(Nx*Nz), 1);

% Relative Permeability

%% Internal Nodes
% Define total head H
Z_array = Z(:);
for i = 1:length(Z_array)
    H(i) = h(Mesh_ref(i)) + Z_array(Mesh_ref(i));
end
%Total Head
Z_array = Z(:);
H = h(Mesh_ref(:)) + Z_array(:);
H_P = H(N_P(:));
H_N = H(N_N(:));
H_E = H(N_E(:));
H_S = H(N_S(:));
H_W = H(N_W(:));

% Calcualting K_tensors on volume faces
Kz_n = (1/2).*(Kzz(N_P(:)) + Kzz(N_P(:)+1));
Kz_s = (1/2).*(Kxx(N_P(:)) + Kxx(N_P(:)-1));
Kx_e = (1/2).*(Kzz(N_P(:)) + Kzz(N_P(:)+Nz));
Kx_w = (1/2).*(Kxx(N_P(:)) + Kxx(N_P(:)-Nz));

[S] = Sat_func(Mesh_ref, ns, h, alpha);
[psi_old] = psi_func(h, psi_res, psi_sat, S, Mesh_ref)
[k] = Relperm_func(ns, ms, h, S, Mesh_ref);

% Calculating k(h) on volume faces
k_P = k(N_P(:));
k_N = k(N_N(:));
k_E = k(N_E(:));
k_S = k(N_S(:));
k_W = k(N_W(:));
avg_type = "upstream";
[k_n, k_s, k_e, k_w] = approximate_k(H_P, H_N, H_E, H_S, H_W, k_P, k_N, k_S, k_E, k_W, avg_type, dz_n, dz_s, dx_e, dx_w);

% Calculating fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using STATIC MESH values for testing purposes%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(Mesh_ref)
    q_e_old(i) = -dx(2,2)/2 * k_e(i) * Kx_e(i) * ((H_E(i) - H_P(i))/10);
    q_n_old(i) = -dz(2,2)/2 * k_n(i) * Kz_n(i) * ((H_N(i) - H_P(i))/5);
    q_w_old(i) = -dx(2,2)/2 * k_w(i) * Kx_w(i) * ((H_P(i) - H_W(i))/10);
    q_s_old(i) = -dz(2,2)/2 * k_s(i) * Kz_s(i) * ((H_P(i) - H_S(i))/5);
end

for i = 1:length(N_P)

    % Define the discretisation for internal nodes    
    F(N_P(i)) = dt/(dx(N_P(i))*dz(N_P(i)))*(theta1*(q_w - q_e) + (1-theta1)*(q_w_old - q_e_old))...
              + dt/(dx(N_P(i))*dz(N_P(i)))*(theta1*(q_s - q_n) + (1-theta1)*(q_s_old - q_n_old))...
              + dt*(theta2*Q(N_P(i)) + (1-theta2)*Q_old(N_P(i)))...
              + psi_old(N_P(i)) - psi(h);
end

%% Bottom Left Corner (0,0) (no flow in -i or -j direction

% Define total head H
    H_P = H(Nz);
    H_N = H(Nz-1);
    H_E = H(2*Nz);
    H_W = 0;
    H_S = 0;
    
    % Determine k_n, k_e
    % k_W = NaN, k_S = NaN
    k_P = k(Nz);
    k_N = k(Nz-1);
    k_E = k(2*Nz);
    k_W = 0;
    k_S = 0;
    avg_type = "arithmetic";
    [k_n, k_e] = approximate_k(H_P, H_N, H_E, H_S, H_W, k_P, k_N, k_S, k_E, k_W, avg_type, dz_n, dz_s, dx_e, dx_w);
    
    % Average hydraulic conductivity of subcontrol volumes
    Kx_e = (1/2).*(Kxx(Mesh_ref(Nz)) + Kxx(Mesh_ref(Nz)+Nz));
    Kx_n = (1/2).*(Kzz(Mesh_ref(Nz)) + Kzz(Mesh_ref(Nz-1)));
        
    % Define flux across control volume faces
    q_e = -(dx_e(1)/2)/2 * k_e * Kx_e * ((H_E - H_P)/dx_e);
    q_n = -(dz_n(1)/2)/2 * k_n * Kz_n * ((H_N - H_P)/dz_n);
    
% Define the discretisation for internal nodes    
    F(Mesh_ref(end, 1)) = dt/(dx(Mesh_ref(end, 1))*dz(Mesh_ref(end, 1)))*(theta1*(-q_e) + (1-theta1)*(-q_e_old))...
              + dt/(dx(Mesh_ref(end, 1))*dz(Mesh_ref(end, 1)))*(theta1*(-q_n) + (1-theta1)*(-q_n_old))...
              + dt*(theta2*Q(Mesh_ref(end, 1)) + (1-theta2)*Q_old(Mesh_ref(end, 1)))...
              + psi_old(Mesh_ref(end, 1)) - psi(Mesh_ref(end, 1));
          
%% Bottom Right Corner (500,0) (no flow in -j direction, CSG dependent in i)

    % Define total head H
    H_P = H(end);
    H_N = H(end - 1);
    H_W = H(end - Nz);
    H_E = 0;
    H_S = 0;
    
    % Determine k_n, k_w
    % k_E = NaN, k_S = NaN
    k_P = k(Mesh_ref(end));
    k_N = k(Mesh_ref(end-1));
    k_W = k(Mesh_ref(end - Nz));
    k_E = 0;
    k_S = 0;
    avg_type = "arithmetic";
    [k_n, k_w] = approximate_k(H_P, H_N, H_E, H_S, H_W, k_P, k_N, k_S, k_E, k_W, avg_type, dz_n, dz_s, dx_e, dx_w);
    
    % Average hydraulic conductivity of subcontrol volumes
    Kx_w = (1/2).*(Kxx(Mesh_ref(end)) + Kxx(Mesh_ref(end-Nz)));
    Kx_n = (1/2).*(Kzz(Mesh_ref(end)) + Kzz(Mesh_ref(end-1)));
        
    % Define flux across control volume faces
    q_w = -(dx_e(1)/2)/2 * k_w * Kx_w * ((H_P - H_W)/dx_w);
    q_n = -(dz_n(1)/2)/2 * k_n * Kz_n * ((H_N - H_P)/dz_n);
    
    if CSG == true
        q_e = Kxx*dH_CSG;
    else
        q_e = 0;
    end

    % Define the discretisation for right boundary nodes  
    F(Mesh_ref(i, end)) = dt/(dx(Mesh_ref(i, end))*dz(Mesh_ref(i, end)))*(theta1*(q_w - q_e) + (1-theta1)*(q_w_old - Kxx*dH_CSG_old))...
          + dt/(dx(Mesh_ref(i, end))*dz(Mesh_ref(i, end)))*(theta1*(q_s - q_n) + (1-theta1)*(q_s_old - q_n_old))...
          + dt*(theta2*Q(Mesh_ref(i, end)) + (1-theta2)*Q_old(Mesh_ref(i, end)))...
          + psi_old(Mesh_ref(i, end)) - psi(Mesh_ref(i, end)); 

%% Top Left Corner (0,100) (river in -i direction, evapotranspiration, rainfall in j)

%% Top Right Corner (500,100) (no flow in i, evapotranspiration, rainfall in j)

%% Left Boundary nodes x = 0 (no flow)

for i = Mesh_ref(2, 1):Mesh_ref(end-1, 1)
    % Define total head H
    H_P = h(Mesh_ref(i, 1)) + Z(i);
    H_N = h(Mesh_ref(i-1, 1)) + Z(i-1);
    H_S = h(Mesh_ref(i+1, 1)) + Z(i+1);
    H_E = h(Mesh_ref(i, 2)) + Z(i);
    
    % Determine k_n, k_s, k_e
    % k_W = NaN
    [k_n, k_s, k_e, ~] = approximate_k(k_P, k_N, k_S, k_E, NaN, avg_type);
    
    % Average hydraulic conductivity of subcontrol volumes
    Kx_e = 1/2*(Kxx(N_P(i) + (Nz)) + Kxx(N_P(i))); % Need to specify Kx_NE, Kx_SE
    Kz_n = 1/2*(Kzz(N_P(i) + 1) + Kzz(N_P(i)));
    Kz_s = 1/2*(Kzz(N_P(i) - 1) + Kzz(N_P(i)));
    
    % Define flux across control volume faces
    q_e = -1/2 * k_e * Kx_e * ((H_E - H_P)/dx_e);
    q_n = -1/2 * k_n * Kz_n * ((H_N - H_P)/dz_n);
    q_s = -1/2 * k_s * Kz_s * ((H_P - H_S)/dx_s);

    % Define the discretisation for left boundary nodes
    if z > 75
        % include river discharge condition
        q_w = -K_R*delta_H_R;
    else
        F(Mesh_ref(i, 1)) = dt/(dx(Mesh_ref(i, 1))*dz(Mesh_ref(i, 1)))*(theta1*(-q_e) + (1-theta1)*(-q_e_old))...
              + dt/(dx(Mesh_ref(i, 1))*dz(Mesh_ref(i, 1)))*(theta1*(q_s - q_n) + (1-theta1)*(q_s_old - q_n_old))...
              + dt*(theta2*Q(Mesh_ref(i, 1)) + (1-theta2)*Q_old(Mesh_ref(i, 1)))...
              + psi_old(Mesh_ref(i, 1)) - psi(Mesh_ref(i, 1));
    end

end

%% Right Boundary nodes x = 500 (no flow)

for i = Mesh_ref(2, end):Mesh_ref(end-1, end)
    
    % Define total head H
    H_P = h(Mesh_ref(i, end)) + Z(i);
    H_N = h(Mesh_ref(i-1, end)) + Z(i-1);
    H_S = h(Mesh_ref(i+1, end)) + Z(i+1);
    H_W = h(Mesh_ref(i, end-1)) + Z(i);
    
    % Determine k_n, k_s, k_w 
    % Note: k_E not defined so set to NaN
    [k_n, k_s, ~, k_w] = approximate_k(k_P, k_N, k_S, NaN, k_W, avg_type);
    
    % Average hydraulic conductivity
    Kx_w = 1/2*(Kxx(N_P(i) - (Nz)) + Kxx(N_P(i)));
    Kz_n = 1/2*(Kzz(N_P(i) + 1) + Kzz(N_P(i)));
    Kz_s = 1/2*(Kzz(N_P(i) - 1) + Kzz(N_P(i)));
    
    % Flux
    q_n = -1/2 * k_n * Kz_n * ((H_N - H_P)/dz_n);
    q_w = -1/2 * k_w * Kx_w * ((H_P - H_W)/dx_w);
    q_s = -1/2 * k_s * Kz_s * ((H_P - H_S)/dz_s);
    
    if Z(i) < 5 && CSG == true
        q_e = Kxx*dH_CSG;
    else
        q_e = 0;
    end

    % Define the discretisation for right boundary nodes  
    F(Mesh_ref(i, end)) = dt/(dx(Mesh_ref(i, end))*dz(Mesh_ref(i, end)))*(theta1*(q_w - q_e) + (1-theta1)*(q_w_old - Kxx*dH_CSG_old))...
          + dt/(dx(Mesh_ref(i, end))*dz(Mesh_ref(i, end)))*(theta1*(q_s - q_n) + (1-theta1)*(q_s_old - q_n_old))...
          + dt*(theta2*Q(Mesh_ref(i, end)) + (1-theta2)*Q_old(Mesh_ref(i, end)))...
          + psi_old(Mesh_ref(i, end)) - psi(Mesh_ref(i, end)); 
end

%%
