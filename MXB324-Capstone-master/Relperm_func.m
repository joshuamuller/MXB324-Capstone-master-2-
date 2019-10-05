function [k] = Relperm_func(ns, ms, h, S)
% Calculates the function k(h), which determines the relative permeability at
% every face.
% Author: Joshua
k = [];
for i = 1:length(ns)
if h(i) < 0
    k(i) = sqrt(S(i)).*(1-(1-S(i)^(1./ms(i))).^(ms(i)))^2;
elseif h(i) >= 0 
    k(i) = 1;
end
end

