function [s, rowsForColumn] = jacobian_directions(N)
% PENDING DOCUMENTATION BY PAUL


Nx = N(1);
Nz = N(2);

area = Nx*Nz;
figure()
hold on
J = diag(ones(area, 1), 0)...
    + diag(ones(area-1, 1), -1)...
    + diag(ones(area-1, 1), 1);

J = J + diag(ones(area - Nx, 1), Nx) + diag(ones(area-Nx, 1), -Nx);

spy(J);

uncaptured = area;
handled = zeros(area, 1);

i = 1;
s = [];

while uncaptured > 0
    vec = zeros(area, 1);
    for col = 1:area
        if handled(col) == 0
            collided = 0;
            for row = 1:area
                if J(row, col) ~= 0
                    for subcol=1:col
                        if col == subcol
                            break;
                        elseif vec(subcol) == 0
                            continue;
                        elseif J(row, subcol) ~= 0
                            collided = 1;
                        end
                        
                    end
                    
                    if vec(row) ~= 0
                        collided = 1;
                    end

                    if collided == 1
                        break; 
                    end
                end
            end
            
            if collided == 0
                for row=1:area
                    if vec(row) == 0
                        continue; 
                    end
                    for j=1:size(s, 2)
                        if s(row, j) ~= 0
                            vec(row) = 0;
                        end
                    end
                end
                
                vec(col) = 1; %vec + J(:, col);
            end
        else
            
        end
    end
    
    for col = 1:area
        if vec(col) ~= 0
            handled(col) = 1; 
            uncaptured = uncaptured - 1;
        end
    end
    
    s(:, i) = vec;
    i = i + 1;
end

rowsForColumn = zeros(area, 5);
for column=1:area
    rows = [5, 2, 3, 1, 1];
    i = 1;
    for row=1:area
        if J(column, row) ~= 0
            rows(1, i) = row;
            i = i + 1;
        end
    end
    
    rowsForColumn(column, :) = rows;
end

end
