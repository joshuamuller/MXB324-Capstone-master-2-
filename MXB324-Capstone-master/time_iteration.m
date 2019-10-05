function [h, yearly] = time_iteration(N, L, Mesh_ref, N_P, N_N, N_S, N_W,...
    N_E, dx, dz, K_xx, K_zz, bore, river, vegetation, rain_percentages, ...
    rainfall, psis, ns, alpha, years)
%SOLVE_SYSTEM Solves for the water content over the given number of years
%
%Outputs: 
%- h = The hydraulic head at each point of the mesh at time zero.
%- yearly = The total water content at the end of each year of simulation.
%
%Inputs:
%- N = The node count on each axis:
%   - (1) Nx
%   - (2) Nz
%- L = Length values of each axis:
%   - (1) Lx
%   - (2) Lz
%- dx = Vector of differences between each node on the x axis
%- dz = Vector of differences between each node on the z axis
%- Ks = All the hydraulic conductivities, as pairs for the x then z axis.
%   - (1) K_allu
%   - (2) K_conf
%   - (3) K_wall
%   - (4) K_volc
%   - (5) K_riv
%- bore = The bore properties:
%   - (1) bore_x
%   - (2) bore_bottom
%   - (3) bore_top
%   - (4) bore_vert_flow_factor = The flow factor of vertical conductivity
%   - (5) bore_rain_percentage = The max percentage of yearly rainfall
%- river = The river properties
%   - (1) riv_b = River bottom (measured from origin)
%   - (2) riv_d = River depth
%   - (3) riv_t = River width (measured from origin)
%- vegetation = The height of each evapotranspiration zone
%   - (1) l_repveg = The repairing vegetation zone (i=1)
%   - (2) l_crop = The crop zone (i=2)
%   - (3) l_plant = The plantation zone (i=3)
%- rain_percentages = The percentage of rain in each zone of the vegetation
%   - (1) rain_repveg
%   - (2) rain_crop
%   - (3) rain_plant
%- rainfall = The rainfall for each day of the simulation
%- psis = The saturated water contents 
%   - (1) psi_allu
%   - (2) psi_conf
%   - (3) psi_wall
%   - (4) psi_volc
%- ns = The n parameter for each area of the aquifer
%   - (1) n_allu
%   - (2) n_conf
%   - (3) n_wall
%   - (4) n_volc
%- alphas = The alpha parameter for each area of the aquifer
%   - (1) alpha_allu
%   - (2) alpha_conf
%   - (3) alpha_wall
%   - (4) alpha_volc
%- years = The number of years the simulation will go over.
%
%Task name: time_iteration
%Authors: Paul
%
%Revisions [Date, people, reason]:
%- 15/09/2019, Paul, Added a lot of the Shamanskii-Newton from Bench3.
%- 21/09/2019, Paul, Reworked for proper system function and Jacobian FD.

ms = 1 - 1./ns;

h = initial_mesh(N, L, bore, river, vegetation, psis);

% TODO initial values for q olds
[q_n_old, q_e_old, q_s_old, q_w_old] = qfunc();

yearly = zeros(length(h), years);

dt = 1;
final_t = years * 365 * 24 * 60 * 60;
t = dt;
seconds_per_year = 365 * 24 * 60 * 60;

[s, rowsForColumn] = jacobian_directions(N);
Jfunc = @(F, Ffunc, h, q_n_old, q_e_old, q_s_old, q_w_old) jacobian_fd(F,...
    Ffunc, h, q_n_old, q_e_old, q_s_old, q_w_old, s, rowsForColumn, N(1));

% This was a random demo function to test my Jacobian re-ordering.
% What's important to observe is that every node is only dependent on
% itself, so the Jacobian it draws is a 1-diagonal matrix. Ours will not
% be like that, ours will have 5 diagonals, but this helped me confirm that
% there's no pollution in the shape of the Jacobian reordering
%Ffunc = @(h) h.*h + 10; 
%x0 = ones(Nx * Nz, 1);
%F = Ffunc(x0);
%J = jacobian_fd(F, Ffunc, x0, s, rowsForColumn, Nx);
%figure()
%spy(J);

year = 0;

while year < years 
    Ffunc = @(h, q_n, q_e, q_s, q_w, S, psi) system_function(h, N, Mesh_ref, ...
        N_P, N_N, N_S, N_W, N_E, dx, dz, L, K_xx, K_zz, bore, river, vegetation, ...
        rain_percentages, rainfall, S, psi_old, psi, q_n_old, q_e_old, q_s_old, ...
        q_w_old, ns, ms, alpha, t, Mesh_ref);
    
    m = 1; % base Newton method for now
    [new_h, q_n, q_e, q_s, q_w, psi, converged] = newton_solve(Ffunc, ...
        Jfunc, h, m, 20, psi_old, q_n_old, q_e_old, q_s_old, q_w_old, alpha);
    
    if converged == false
        t = t - dt;
        dt = dt/2;
    else
        h = new_h;
        dt = dt * 3;
        q_n_old = q_n;
        q_e_old = q_e;
        q_s_old = q_s;
        q_w_old = q_w;
        psi_old = psi;
    end
    
    
    % Apply time step
    t = t + dt;
    
    new_year = mod(t, seconds_per_year);
    if new_year > year
        year = new_year; 
        yearly(:, year) = h;
    end
    
end

end