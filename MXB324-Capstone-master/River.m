function [D_H] = River(h, Z, K_R, dz)
% Calculates river discharge at the top left hand corner of the domain

H_R = 85; D_xR = 50; R_t = 50; R_b = 75; R_d = 25; 

for i = 1:(25/dz(2,2)+1)
if h(i) + Z(i) <= H_R
    D_H(i) = (H_R - (h(i)+Z(i)))/D_xR;
elseif h(i) + Z(i) <= R_d
    D_H(i) = (H_R - R_b)/D_xR;
end
end
end

