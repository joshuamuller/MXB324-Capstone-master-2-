function [psi_old] = psi_func(h, psi_res, psi_sat, S)
% Calculates the function Psi(h), which determines the water content at
% every node.
% Author: Joshua
psi_old = [];
for i = 1:length(h)
    if h(i) < 0
    	psi_old(i) = psi_res(i) + S(i).*(psi_sat(i) - psi_res(i));
    elseif h(i) >= 0
        psi_old(i) = psi_sat(i);
    end
end
end
