function [h] = hScoreCalc3D(n, Target, Risk, Vox)
% This function calculates the hscore for a particular node (n = [x,y,z]), 
% this is a heuristic function which estimates the cost to the target.
% Inputs -  n - node [x,y,z]km
%           Target - Target coordinates [x,y,z]km
%           Risk  - Envionment risk (3D matrix)
%           Vox - Voxel size [x,y,z]km

%Nodal distance from node to Target
dx = (Target(1)-n(1));
dy = (Target(2)-n(2));
dz = (Target(3)-n(3));
d = (dx^2 + dy^2 + dz^2)^0.5;
if d > 0
    %Cumulative Density to Target
    N = ceil(d/(2^0.5)); %Points interpolated at
    for i = 1:N
        x = n(1) + (i-0.5)*dx/N;
        y = n(2) + (i-0.5)*dy/N;
        z = n(3) + (i-0.5)*dz/N;
        density(i) = interp3(Risk,x,y,z); %Density at each point
    end
    
    %Actual distance
    Dx = (Target(1)-n(1))*Vox(1);
    Dy = (Target(2)-n(2))*Vox(2);
    Dz = (Target(3)-n(3))*Vox(3);
    D = (Dx^2 + Dy^2 + Dz^2)^0.5;
    CD = mean(density)*D;
    h = sum(CD); %Cumulative density to target
else
    h = 0;
end