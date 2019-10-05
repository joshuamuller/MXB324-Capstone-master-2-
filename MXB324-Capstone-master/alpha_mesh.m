function [alpha] = alpha_mesh(alpha_val, X, Z, Nx, Nz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:Nz
    for j = 1:Nx
        if X(1,j) <= 50 && Z(i,1) <= 100 && Z(i,1) >=40 %Alluvium
            alpha(i,j) = alpha_val(1);      
        elseif X(1,j) >= 50 && X(1,j) <= 300 && Z(i,1) <= 100 && Z(i,1) >=50
            alpha(i,j) = alpha_val(1);        
        end
        if X(1,j) >= 300 && X(1,j) <=500 && Z(i,1) <= 100 && Z(i,1) >=50 %Main Range Volcanics
            alpha(i,j) = alpha_val(2);        
        end
        if X(1,j) <= 500 && Z(i,1) <= 100 && Z(i,1) <=40 %Walloon Coal Measures
            alpha(i,j) = alpha_val(3);
            
        elseif X(1,j) <= 500 && X(1,j)>= 300 && Z(i,1) <= 50 && Z(i,1) >=40
            alpha(i,j) = alpha_val(3);         
        end
        if X(1,j) >= 50 && X(1,j) <= 350 && Z(i,1) >=40 && Z(i,1) <= 50 %Confining Layer
            alpha(i,j) = alpha_val(4);
        end
    end
end
alpha = alpha(:);
end