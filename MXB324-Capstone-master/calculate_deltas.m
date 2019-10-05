function [dx, dz] = calculate_deltas(method, Nx, Nz, Lx, Lz, bore, ...
                                                        refiner_weights)
%CALCULATE_DELTAS Creates the mesh deltas
%
%Outputs:
%- dx = The distance between two nodes along the x axis, a function of x.
%- dz = The distance between two nodes along the z axis, a function of z.
%
%Inputs:
%- method = The method to use ('uniform' or 'refined')
%- Nx = The number of nodes on each row
%- Nz = The number of nodes on each column
%- Lx = Horizontal length of the problem
%- Lz = Vertical height of the problem
%- bore = Container of [bore_x, bore_bottom, bore_top].
%   - (1) bore_x = The x coordinate of the bore
%   - (2) bore_bottom = The lowest z coordinate of the bore
%   - (3) bore_top = The highest z coordinate of the bore
%
%The mesh is a one dimensional array but the geometry of the mesh is 2D 
%with, potentially, varying distances between each node point for 
%vertical and horizontal separately. This method takes the method 
%(either uniform or refined), and returns the array of distances for both 
%the x and z axes.
%
%Task name: mesh
%Authors: Paul
%
%Revisions [Date, people, reason]:
%- 15/09/2019, Paul, Added a lot of the Shamanskii-Newton from Bench3.

bore_x = bore(1);
bore_bottom = bore(2);
bore_top = bore(3);

if strcmp(method, 'uniform')
    dx = ones(Nx+1, 1) .* (Lx/(Nx+1));  
    dz = ones(Nz+1, 1) .* (Lz/(Nz+1));  
elseif strcmp(method, 'refined')
    bore_x_rel = bore_x/Lx;
    bore_min_z_rel = bore_bottom / Lz;
    bore_max_z_rel = bore_top / Lz;
    
    % In order: first 5% of nodes on the x axis, 5% before and after bore,
    %           last 5% of nodes on the x axis.
    tight_areas_x = [1, floor(Nx/15), 0 ;
                     floor([0.8, 1.2] .* bore_x_rel .* Nx), 1];
    
    dx = refiner(Nx, Lx, refiner_weights, tight_areas_x);
    
    x_sum = 0;
    x_transitions = [50, 300, 350];
    for i=1:size(dx, 1)
        for j=1:length(x_transitions)
            tr = x_transitions(j); 
            if x_sum < tr && x_sum + dx(i) > tr
                carryOverSize = x_sum + dx(i) - tr;
                if carryOverSize > dx(i)/2
                    carryOverSize = dx(i) - carryOverSize;
                    dx(i-1) = dx(i-1) + carryOverSize;
                    dx(i) = dx(i) - carryOverSize;
                    x_sum = x_sum + carryOverSize;
                else
                    dx(i) = dx(i) - carryOverSize;
                    dx(i+1) = dx(i+1) + carryOverSize;
                end
            end
        end
        x_sum = x_sum + dx(i);
    end
    
    % In order: last 5% of nodes on the z axis, 5% of nodes prior to and
    %           after the bore range of the z axis
    tight_areas_z = [floor(0.9 * Nz), Nz, 2 ;
                     floor(0.9 * bore_min_z_rel * Nz), ceil(1.1 * bore_max_z_rel * Nz), 1];
    dz = refiner(Nz, Lz, refiner_weights, tight_areas_z);
    
    z_sum = 0;
    z_transitions = [40, 50];
    for i=1:size(dz, 1)
        for j=1:length(z_transitions)
            tr = z_transitions(j); 
            if z_sum < tr && z_sum + dz(i) > tr
                carryOverSize = z_sum + dz(i) - tr;
                if carryOverSize > dz(i)/2
                    carryOverSize = dz(i) - carryOverSize;
                    dz(i-1) = dz(i-1) + carryOverSize;
                    dz(i) = dz(i) - carryOverSize;
                    z_sum = z_sum + carryOverSize;
                else
                    dz(i) = dz(i) - carryOverSize;
                    dz(i+1) = dz(i+1) + carryOverSize;
                end

                break;
            end
        end
        z_sum = z_sum + dz(i);
    end
end

% Comment next bit out when you don't want to plot
% Pssst, Ctrl-R to comment, Ctrl-T to uncomment.
figure();
xlim([0, Lx]);
ylim([0, Lz]);
hold on;
   
x_sum = 0;
z_sum = 0;
   
plot([0, 0], [0, Lz], 'b');
for i=1:Nx
   x_sum = x_sum + dx(i);
   plot([x_sum, x_sum], [0, Lz], 'b');
end
   
plot([0, Lx], [0, 0], 'r');
for i=1:Nz
   z_sum = z_sum + dz(i);
   plot([0, Lx], [z_sum, z_sum], 'r');
end

x = zeros(Nx, 1);
sum = 0;
for xi = 1:Nx
    sum = sum + dx(xi);
    x(xi) = sum;
end

z = zeros(Nz, 1);
sum = 0;
for zi = 1:Nz
    sum = sum + dz(zi);
    z(zi) = sum;
end


figure()
hold on;
for xi=1:Nx
    scatter(x(xi) * ones(size(z)), z, 8, 'b', 'filled');
    
    for zi = 1:Nz
        minx = x(xi) - dx(xi)/2;
        maxx = x(xi) + dx(xi+1)/2;
        minz = z(zi) - dz(zi)/2;
        maxz = z(zi) + dz(zi+1)/2;
        
        plot([minx, minx, maxx, maxx, minx],  ...
             [minz, maxz, maxz, minz, minz], 'r--');
    end
end

plot([50, 50], [0, 100], 'k-');
plot([300, 300], [0, 100], 'k-');
plot([350, 350], [0, 100], 'k-');

plot([0, 500], [40, 40], 'k-');
plot([0, 500], [50, 50], 'k-');

end
