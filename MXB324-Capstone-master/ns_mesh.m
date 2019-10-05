function [ns] = ns_mesh(ns_val, X, Z, Nx, Nz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:Nz
    for j = 1:Nx
        if X(1,j) <= 50 && Z(i,1) <= 100 && Z(i,1) >=40 %Alluvium
            ns(i,j) = ns_val(1);      
        elseif X(1,j) >= 50 && X(1,j) <= 300 && Z(i,1) <= 100 && Z(i,1) >=50
            ns(i,j) = ns_val(1);        
        end
        if X(1,j) >= 300 && X(1,j) <=500 && Z(i,1) <= 100 && Z(i,1) >=50 %Main Range Volcanics
            ns(i,j) = ns_val(2);        
        end
        if X(1,j) <= 500 && Z(i,1) <= 100 && Z(i,1) <=40 %Walloon Coal Measures
            ns(i,j) = ns_val(3);
            
        elseif X(1,j) <= 500 && X(1,j)>= 300 && Z(i,1) <= 50 && Z(i,1) >=40
            ns(i,j) = ns_val(3);         
        end
        if X(1,j) >= 50 && X(1,j) <= 350 && Z(i,1) >=40 && Z(i,1) <= 50 %Confining Layer
            ns(i,j) = ns_val(4);
        end
    end
end
ns = ns(:);
end