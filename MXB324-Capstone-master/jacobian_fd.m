function J = jacobian_fd(F, Ffunc, h, q_n_old, q_e_old, q_s_old, q_w_old,...
    s, rowsForColumn, Nx)
%PENDING DOCUMENTATION BY PAUL

N = size(F, 1);
J = zeros(N, N);
sqeps = sqrt(eps);

diff = norm(h, 2);
if diff == 0
    diff = diff * sqeps;
else
    diff = sqeps;
end

Jnew = zeros(N, size(s, 2));
for j = 1:size(s, 2)
    e = s(:, j);
    enorm = norm(e, 2);
    diff = diff / enorm;
    Jnew(:, j) = (Ffunc(h + diff * e, q_n_old, q_e_old, q_s_old, q_w_old) - F)/diff;
end

for column=1:N
    % Find which unit vector took this column
    for i=1:size(s, 2)
        if s(column, i) ~= 0
            rs = rowsForColumn(column, :);
            
            if column == 1 || column == N
                rs = rs([1, 2, 3]);
            elseif column <= Nx || column > N - Nx
                rs = rs([1, 2, 3, 4]);
            end
            
            % Oh good, now use that column of Jnew since that's where it is
            J(rs, column) = Jnew(rs, i); 
            break;
        end
    end
end





end