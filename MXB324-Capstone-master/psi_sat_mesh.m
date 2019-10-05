function [psi_sat] = psi_sat_mesh(psi_sat_val, X, Z, Nx, Nz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:Nz
    for j = 1:Nx
        if X(1,j) <= 50 && Z(i,1) <= 100 && Z(i,1) >=40 %Alluvium
            psi_sat(i,j) = psi_sat_val(1);      
        elseif X(1,j) >= 50 && X(1,j) <= 300 && Z(i,1) <= 100 && Z(i,1) >=50
            psi_sat(i,j) = psi_sat_val(1);        
        end
        if X(1,j) >= 300 && X(1,j) <=500 && Z(i,1) <= 100 && Z(i,1) >=50 %Main Range Volcanics
            psi_sat(i,j) = psi_sat_val(2);        
        end
        if X(1,j) <= 500 && Z(i,1) <= 100 && Z(i,1) <=40 %Walloon Coal Measures
            psi_sat(i,j) = psi_sat_val(3);
            
        elseif X(1,j) <= 500 && X(1,j)>= 300 && Z(i,1) <= 50 && Z(i,1) >=40
            psi_sat(i,j) = psi_sat_val(3);         
        end
        if X(1,j) >= 50 && X(1,j) <= 350 && Z(i,1) >=40 && Z(i,1) <= 50 %Confining Layer
            psi_sat(i,j) = psi_sat_val(4);
        end
    end
end
psi_sat = psi_sat(:);
end