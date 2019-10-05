function [h, q_n, q_e, q_s, q_w, psi, succeeded] = newton_solve(Ffunc, Jfunc, ...
    h_old, m, max_iters, psi_old, q_n_old, q_e_old, q_s_old, q_w_old, alpha)
%NEWTON_SOLVE Performs Shamanskii-Newton iteration to approximate 
% the next h. It can fail, in which case succeeded will be false.
%
%Outputs:
%- h = The next step of h.
%- succeeded = true if it converged in time, false otherwise.
%
%Inputs:
%- Ffunc = The function to approximate the zero of, takes h, q_n_old etc, S, psi_old, psi.
%- Jfunc = The function for calculating the Jacobian, takes F, Ffunc, h
%- h_old = The previous solution
%- deltas = The delta values for the x and y axis
%   - (1) = The vector containing all delta values for the x axis
%   - (2) = The vector containing all delta values for the z axis
%- m = The Jacobian refresh interval. As this increases, the performance
%  increases at the cost of accuracy. m=1 is the regular Newton method.
%- max_iters = The maximum number of iterations it is allowed.
%
%Task name: newton
%Authors: Paul
%
%Revisions [Date, people, reason]:
%- 15/09/2019, Paul, Added a lot of the Shamanskii-Newton from Bench3.
%- 21/09/2019, Paul, Reworked for proper system function and Jacobian FD.

% Set relative and Absolute Tolerances

S = sat_func(Mesh_ref, ns, h_old, alpha);
h = h_old; % Our initial newton guess is the previous h value

psi = psi_old; % Our initial guess is the previous h value, so psis will be the same for the first guess


F = Ffunc(h_old, q_n_old, q_e_old, q_s_old, q_w_old, S, psi);

%(h, q_n, q_e, q_s, q_w, S, psi_old, psi
k = 1; 
err = Inf; 
errold=norm(F,2);
errvecS = zeros(max_iters,1);
errvecS(1) = errold;
tola = 1e-10;
tolr = 1e-12;
tol = tola + tolr*errvecS(1); % Scale Relative tolerance by initial error.


% Newton Iteration
while err > tol && k < max_iters
    if mod(k, m) ==0
        k = 0; 
        J = Jfunc(F, Ffunc, h, q_n_old, q_e_old, q_s_old, q_w_old); 
    end
    
    dh = J\(-F);
    h = h + dh;
    
    % Update S and psi, q will probably end up out here as well
    S = sat_func(Mesh_ref, ns, h, alpha);
    psi = psi_func(h, psi_res, psi_sat, S, Mesh_ref);
    
    
    % TODO
    % Put [q_n, q_e, q_s, q_w] = q_func(...)
    
    % Calculate how wrong this guess was
    [F, q_n, q_e, q_s, q_w] = Ffunc(h, q_n_old, q_e_old, q_s_old, q_w_old, S, psi);
    err = norm(F,2); 
    k = k + 1;
end

succeeded = err > tol;

end

