function [CD_actual,Distance] = ActualCD3D(Vox,Size, Map,Path)
% This function calculates the actual cumulative desity of a path.
% Inputs -  Vox - Voxel size [x,y,z]km
%           Start - Start coordinates [x,y,z]km
%           Target - Target coordinates [x,y,z]km
%           Map  - Envionment risk (3D matrix)
%           Path - Coordinates of path [x,y,z]km

MapRes = size(Map)-1;      %Resolution of Map
S_Map = MapRes./Size;      %Map scale

Path_Map = Path.*S_Map + 1;%Path scaled to Map size

CD_actual = 0;             %Cumulative densitiy at start of path 
Distance = 0;

%Number of interpolated points based on resolution in x
P = S_Map(1)*Vox(1);
pi = Vox.*MapRes./(Size*P);
%loop through path
for n = 2:size(Path,1)
    %Step distance
    dd = Path(n,:)-Path(n-1,:);
    d = (sum(dd.^2))^0.5;
    %Scaled distance on map based on resolution
    dd_Map = dd./Vox;
    %Interpolate points between Step

    for i = 1:P
        p = Path_Map(n-1,:) + (i-0.5)*pi.*dd_Map;
        Risk(i) = interp3(Map,p(2),p(1),p(3));
    end
    Risk = mean(Risk); %Mean Risk
    CD_actual = CD_actual + Risk*d;
    Distance = Distance + d;
end

 