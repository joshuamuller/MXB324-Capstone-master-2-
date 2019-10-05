function [dx, dz, dx_node, dz_node] = Deltas(Nx,Nz, X, Z)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Calculate distance between nodes

for i = 1:length(X(1,:))-1
    for j = 1:length(Z(:,1))-1
        dx_node(i,j) = abs(X(j,i) - X(j,i+1));
        dz_node(i,j) = abs(Z(j,i) - Z(j+1,i));
    end
end

% Calculate control volume face lengths

for j = 1:Nz-1
    for i = 2:Nx-1
        dx(j,i) = abs(X(j,i-1) - X(j,i))/2 + abs(X(j,i) - X(j,i+1))/2;
    end
end
dx(:,1) = abs(X(1,1)-X(1,2))/2; dx(:,Nx) = abs(X(1,Nx-1)-X(1,Nx))/2;
for j = 1:Nx-1
    for i = 2:Nz-1
        dz(i,j) = abs(Z(i-1,j) - Z(i,j))/2 + abs(Z(i+1,j) - Z(i,j))/2;
    end
end
dz(1,:) = abs(Z(1,1)-Z(2,1))/2; dz(Nz,:) = abs(Z(Nz,1)-Z(Nz-1,1))/2;
end

