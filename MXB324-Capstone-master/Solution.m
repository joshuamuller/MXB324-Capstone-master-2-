%% Solution file
%Task name: variables
%Authors: Paul
%Revisions: [Date, people, reason];
%- 21/09/2019, Paul, Reworked for proper system function and Jacobian FD.

%% Clean workspace

format compact;
clear;
clc;
close all;

%% Parameters

% Behaviour
rainfall_model = 'predicted';   % or 'drought'/'flood'
years = 50;                     % how many years in the simulation
drought_years = [];             % year numbers for drought periods
flood_years = [];               % year numbers for flood periods

% Indexing Conveniences
ALLUVIUM = 1;
CONFINEMENT = 2;
WALLOON = 3;
VOLCANICS = 4;

% n's
ns_val = [1.51 ; 1.3954 ; 2.239 ; 2.0];

% alpha's
alpha_val = [1.43 ; 1.04 ; 2.8 ; 2.5];

% Hydraulic Conductivity (Permeability)
K_xx = [2.6 0.08 3.9 4.3];
K_zz = [0.91 0.016 1.17 0.21];

% Water Contents
psi_res_val = [0.01 ; 0.106 ; 0.01 ; 0.01];
psi_sat_val = [0.33 ; 0.4686 ; 0.1 ; 0.075];

% Dimensions (origin is bottom left corner)
Lx = 500;                   % The width of the aquifer (x)
Lz = 100;                   % The height of the aquifer (z)
L = [Lx, Lz];

% River Properties
riv_b = 75;                 % River bottom (measured from origin)
riv_d = 25;                 % River depth
riv_t = 50;                 % River width (measured from origin)   
river_properties = [riv_b, riv_d, riv_t];

% Height of each evapotranspiration zone
l_repveg = 15;              % Repairing vegetation zone (i=1)
l_crop = 5;                 % Crop zone (i=2)
l_plant = 10;               % Plantation zone (i=3)
vegetation = [l_repveg, l_crop, l_plant];

% Percentage of rainfall at the surface
rain_repveg = 0.06;
rain_crop = 0.035;
rain_plant = 0.045;
rain_percentages = [rain_repveg, rain_crop, rain_plant];

% Bore Properties
bore_x = 100;               % Horizontal coordinate of the bore
bore_bottom = 55;           % Vertical coordinate of the bottom of the bore
bore_top = 75;              % Vertical coordinate of the top of the bore
bore_vert_flow_factor = 75; % Factor for vertical hydraulic conductivity
bore_rain_percentage = 0.5; % The max percentage of yearly rainfall
bore = [bore_x, bore_bottom, bore_top, bore_vert_flow_factor, ...
    bore_rain_percentage];

% Mesh Control
Nx = 51;                   % The number of horizontal mesh points
Nz = 21;                   % The number of vertical mesh points

[X, Z, Mesh_ref, N_P, N_N, N_S, N_W, N_E] = Mesh_creator(Nx,Nz);
[dx, dz, dx_node, dz_node] = Deltas(Nx,Nz, X, Z);
[Kxx, Kzz] = K_Tensor(K_xx,K_zz, X, Z, Nx, Nz, dx, dz);
[alpha] = alpha_mesh(alpha_val, X, Z, Nx, Nz);
[psi_res] = psi_res_mesh(psi_res_val, X, Z, Nx, Nz);
[psi_sat] = psi_sat_mesh(psi_sat_val, X, Z, Nx, Nz);
[ns] = ns_mesh(ns_val, X, Z, Nx, Nz);
ms = 1-1./ns;


% Calculating initial head profile for testing purposes
htop = -10;
hbot = -5;
L2 = 100;
h = zeros(Nz, Nx);
for i = 1:Nz
    for j = 1:Nx
        h(i,j) = hbot + (htop - hbot)*Z(i,j)/L2;
    end
end
h = h(:);

% Define total head H
Z_array = Z(:);
for i = 1:length(Z_array)
    H(i) = h(Mesh_ref(i)) + Z_array(Mesh_ref(i));
end

%% Initial Fluxes

%Total Head
Z_array = Z(:);
H = h(Mesh_ref(:)) + Z_array(:);
[S] = Sat_func(ns, h, alpha);
[psi_old] = psi_func(h, psi_res, psi_sat, S)
[k] = Relperm_func(ns, ms, h, S);

[Kxx_n, Kxx_e, Kxx_s, Kxx_w, Kzz_n, Kzz_e, Kzz_s, Kzz_w] = K_Tensor(K_xx,K_zz, X, Z, Nx, Nz, dx, dz)
avg_type = "harmonic"
[k_n, k_s, k_e, k_w] = approximate_k(H, k, avg_type, dz, dx, Nz, Nx);

%%

% The idea here is that the first number divided by the second number
% is the ratio of refined deltas to unrefined deltas.
refiner_weights = [1, 3];

% [dx, dz] = calculate_deltas('refined', Nx, Nz, Lx, Lz, ...
%                             [bore_x, bore_bottom, bore_top], ...
%                             refiner_weights);

% Rainfall
rainfall = predict_rainfall(years, rainfall_model, drought_years, flood_years);

%% Begin model
                
[h, totals] = time_iteration([Nx, Nz], L, X, Z, Mesh_ref, N_P, N_N, N_S, N_W, N_E, dx, dz, Kxx, Kzz, bore, river_properties, vegetation,...
                 rain_percentages, rainfall, [psi_res, psi_sat], ns, alpha, years, Mesh_ref);
% 
% plot_head(h, dx, dz, years);
% plot_totals(totals, years);
