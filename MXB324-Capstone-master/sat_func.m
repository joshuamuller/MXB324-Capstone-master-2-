function [S] = Sat_func(ns, h, alpha)
% Calculates the function S(h), which determines the relative saturation at
% every node.
% Author: Joshua
ms = 1 - 1./ns;
S = [];
for i = 1:length(ns)
if h(i) < 0
    S(i) = ((1 + (-alpha(i).*h(i)).^(ns(i)))).^(-ms(i));
elseif h(i) >= 0 
    S(i) = 1;
end
end
end

