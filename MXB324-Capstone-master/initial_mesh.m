function h = initial_mesh(N, L, bore, river, vegetation, psis)
%INITIAL_MESH Gets the initial values for the water table. This is t=0.
%
%Outputs: 
%- h = The water head at each point of the mesh at time zero.
%
%Inputs:
%- N = The node count on each axis:
%   - (1) Nx
%   - (2) Nz
%- L = Length values of each axis:
%   - (1) Lx
%   - (2) Lz
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
%- psis = The saturated water contents 
%   - (1) psi_allu
%   - (2) psi_conf
%   - (3) psi_wall
%   - (4) psi_volc
%
%Task name: initial_mesh
%Authors: Paul
%
%Revisions [Date, people, reason]:
%-

h = ones(N(1) * N(2), 1);

end

