function [g] = gScoreCalc3D(n, Current, gScore, Risk, Vox)
% This function calculates the gscore for a particular node (n = [x,y,z]), 
% this is the cost of the start of the path to the node
% Inputs -  n - next node [x,y,z]km
%           Current - Current node [x,y,z]km
%           gScore - Matrix of gSCore Values (3D matrix)
%           Risk  - Envionment risk (3D matrix)
%           Vox - Voxel size [x,y,z]km

%Distance from position to node
dx = (n(1)-Current(1))*Vox(1);
dy = (n(2)-Current(2))*Vox(2);
dz = (n(3)-Current(3))*Vox(3);
d = (dx^2 + dy^2 + dz^2)^0.5;

%Average Density across voexl
x = (n(1)+Current(1))/2; %midpoint x coordinate
y = (n(2)+Current(2))/2; %midpoint y coordinate
z = (n(3)+Current(3))/2; %midpoint z coordinate

density = interp3(Risk,x,y,z);
%Cumulative Density from position to node
CD = d * density;
% g Score
g = CD + gScore(Current(2),Current(1),Current(3));
end

