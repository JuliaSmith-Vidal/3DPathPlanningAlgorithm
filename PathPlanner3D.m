function [Results] = PathPlanner3D(Vox, Start, Target, Size, Ground_Risk,Air_Risk, hw)
%PathPlanner - Function to determin optimum path through a 3D environment
%              based on minimising the cumulative risk.

% author : Julia Smith-Vidal
% email : jksv@hotmail.co.uk
% Inputs -  Vox - Voxel size [x,y,z]km
%           Start - Start coordinates [x,y,z]km
%           Target - Target coordinates [x,y,z]km
%           Size  - Envionment dimensions [x,y,z]km
%           Ground_Risk - Risk Data for environment (3D matrix)
%           Air_Risk - Risk Data for environment (3D matrix)
%           hw - heuristic weiting (between 0 and 1 recomended)
% Output    Results - Cumulative density, Distance, Error, Time

tstart = tic; %start time
Map = Ground_Risk + Air_Risk;

nodes = round(Size./Vox + 1);     %Number of nodes    
nodes = [nodes(2), nodes(1), nodes(3)];  %[row,col,pages][y,x,z]

%Check resolution of Map not exeeded
if nodes(1) > size(Map,1)
    'Resolution greater than map resolution'
end
%Check resolution gives an integer of nodes
if any(mod(Size./Vox,1) >0)
    'Resolution Changed' 
end

%%Resize enviornment
Vox = Size./([nodes(2), nodes(1), nodes(3)]-1); %Update resolution based on number of nodes
n_Start = round(Start./Vox + 1);      %Start node
n_Target = round(Target./Vox + 1);    %Target node

MapRes = size(Map)-1;              %Resolution of Map (y,x,z)
S = MapRes./(nodes-1);             %Scale of conversion

x = [1:S(2):size(Map,2)]; %x node positions 
y = [1:S(1):size(Map,1)]; %y node positions
z = [1:S(3):size(Map,3)]; %z node positions
[X,Y,Z] = meshgrid(x,y,z);
Risk = interp3(Map,X,Y,Z);     %Interpolate Map to find risk at resolution


%%Initial set up
gScore = zeros(nodes);
hScore = zeros(nodes);
fScore_Open = zeros(nodes); %Set of nodes to be evaluated
fScore_Closed = zeros(nodes); %Set of nodes already evaluated
Parent = cell(nodes);

%%Initial Values
h = hScoreCalc3D(n_Start, n_Target, Risk, Vox); %Initial heuristic estimation
hScore(n_Start(2),n_Start(1),n_Start(3)) = h;
fScore_Open(n_Start(2),n_Start(1),n_Start(3)) = h;

%%Loop until target is reached
while fScore_Closed(n_Target(2), n_Target(1), n_Target(3)) == 0
    %Determin node with lowest f score
    [fminy,fminx,fminz] = ind2sub(size(fScore_Open),find(fScore_Open == min(fScore_Open(fScore_Open>0))));
    for p = 1 : length(fminx)
        Current = [fminx(p),fminy(p),fminz(p)];
        %Set Current to Closed
        fScore_Closed(Current(2),Current(1),Current(3)) = fScore_Open(Current(2),Current(1),Current(3));
        fScore_Open(Current(2),Current(1),Current(3)) = 0;
        
        % For the start node look at all surronding nodes
        if isequal(Current,n_Start)
            for x = [-1:1:1]
                for y = [-1:1:1]
                    for z = [-1:1:1]
                        %Check path does not go vertically up or down
                        if x == 0 && y == 0 && z ~= 0
                        else
                            %Position to be investigated
                            n = Current + [x , y, z];
                            %Check position is on grid, is a closed node
                            if all(n > 0) && all(n < [nodes(2),nodes(1),nodes(3)]+1) && fScore_Closed(n(2),n(1),n(3)) == 0
                                %Calculate g Score
                                g = gScoreCalc3D(n, Current, gScore, Risk, Vox);
                                %Check if position has a g value or if new value is less
                                if gScore(n(2),n(1),n(3)) == 0 || gScore(n(2),n(1),n(3)) > g
                                    %Assign g Score
                                    gScore(n(2),n(1),n(3)) = g;
                                    if hScore(n(2),n(1),n(3)) == 0
                                        %Assign h Score
                                        h = hScoreCalc3D(n, n_Target, Risk, Vox);
                                        hScore(n(2),n(1),n(3)) = h;
                                    end
                                    %Calculate f Score
                                    fScore_Open(n(2),n(1),n(3)) = g + h*hw;
                                    Parent(n(1),n(2),n(3)) = {Current};
                                end
                            end
                        end
                    end
                end
            end      
   
        else
            %For nodes afte the start calculate current angles of flight
            P = Parent{Current(1),Current(2),Current(3)};
            dx = Current(1) - P(1);
            dy = Current(2) - P(2);
            dz = Current(3) - P(3);
            az = atan2(dy,dx); %Azemuth
            el = atan2(dz,((dx^2+dy^2)^0.5)); %Elevation
            %Investigate nodes within a 45 degree turn angle
            for azc = [-pi/4,0,pi/4]
                for z = [-1,0,1]
                    az_new = az + azc; %New azemuth = Old azemuth + change
                    if mod(az_new,90) == 0
                        d = 1;
                    else
                        d = 2^0.5;
                    end
                    x = round(d*cos(az_new));
                    y = round(d*sin(az_new));

                    %Position to be investigated
                    n = Current + [x , y, z];
                    %Check position is on grid, is a closed node
                    if all(n > 0) && all(n < [nodes(2),nodes(1),nodes(3)]+1) && fScore_Closed(n(2),n(1),n(3)) == 0
                        %Calculate g Score
                        g = gScoreCalc3D(n, Current, gScore, Risk, Vox);
                        %Check if position has a g value or if new value is less
                        if gScore(n(2),n(1),n(3)) == 0 || gScore(n(2),n(1),n(3)) > g
                            %Assign g Score
                            gScore(n(2),n(1),n(3)) = g;
                            if hScore(n(2),n(1),n(3)) == 0
                                %Assign h Score
                                h = hScoreCalc3D(n, n_Target, Risk, Vox);
                                hScore(n(2),n(1),n(3)) = h;
                            end

                            %Calculate f Score
                            fScore_Open(n(2),n(1),n(3)) = g + h*hw;
                            Parent(n(1),n(2),n(3)) = {Current};
                        end
                    end
                end
            end
        end
    end                
end

%%Results
%Reconstruct node Path
n_Path = [n_Target(1),n_Target(2),n_Target(3)];
s = 1; %Steps
%loop until Path meets Start
while ~ismember(n_Start,n_Path,'rows')
    n_Path = [n_Path;Parent{n_Path(s,1),n_Path(s,2),n_Path(s,3)}];
    s = s + 1;
end

time = toc(tstart); %Time taken

CD = gScore(n_Target(2),n_Target(1),n_Target(3)); %Cumulative density

%Scale results for mesh size
Path = flip((n_Path-1)).*Vox; %Path scaled [x, y, z]

%Recalcuate CD
[CD_actual,Distance] = ActualCD3D(Vox,Size, Map,Path);
Error = abs(CD - CD_actual)/CD_actual;
    

Results = [CD_actual,Distance, Error, time];

%%Plot Results
%Calculate azimuth and Elevation
dx = Path(2:end,1) - Path(1:end-1,1); %dx of each step
dy = Path(2:end,2) - Path(1:end-1,2); %dy of each step
dz = Path(2:end,3) - Path(1:end-1,3); %dz of each step
gd = cumsum((dx.^2 + dy.^2).^0.5);    %Ground Distance
az = atan2d(dy,dx);                   %Azemuth (degrees)
el = atan2d(dz,((dx.^2+dy.^2).^0.5)); %Elevation (degrees)
azc = az(2:end)- az(1:end-1);         %Azemuth variation
elc = el(2:end)- el(1:end-1);         %Elevation variation


%%Iso Plot 3D 
figure
subplot(3,3,[1 2 4 5 7 8])
%Plot environment
mshx = [0:Size(1)/MapRes(2):Size(1)];
mshy = [0:Size(2)/MapRes(1):Size(2)];
mshz = [0:Size(3)/MapRes(3):Size(3)];
% Air Riks
p = patch(isosurface(mshx,mshy,mshz(1:end),Air_Risk(:,:,1:end),0.005));
isonormals(mshx,mshy,mshz,Air_Risk,p)
p.FaceColor = 'blue';
p.EdgeColor = 'none';
p.FaceAlpha = 0.2;
hold on 
% Ground Riks
p = patch(isosurface(mshx,mshy,mshz(4:end),Ground_Risk(:,:,4:end),0.01));
isonormals(mshx,mshy,mshz,Ground_Risk,p)
p.FaceColor = 'green';
p.EdgeColor = 'none';
p.FaceAlpha = 0.2;
p = patch(isosurface(mshx,mshy,mshz(1:end),Ground_Risk(:,:,1:end),10));
isonormals(mshx,mshy,mshz,Ground_Risk,p)
p.FaceColor = 	[0,0.5,0];
p.EdgeColor = 'none';
p.FaceAlpha = 0.4;

view([-5 -4 2.5]);
camlight 
lighting gouraud
grid on
hold on
ylabel('y (km)','fontsize',20)
xlabel('x (km)','fontsize',20)
zlabel('z (km)','fontsize',20)
axis([0 8 0 13 0 1])
daspect([1 1 0.5])
ax = gca;
ax.FontSize = 16; 
%Plot Route
plot3(Path(:,1),Path(:,2),Path(:,3),'r.-','LineWidth',1.5,'MarkerSize',13)

x = [0:size(Map,2)-1]*Size(1)/(size(Map,2)-1); %x node positions 
y = [0:size(Map,1)-1]*Size(2)/(size(Map,1)-1); %y node positions
z = [0:size(Map,3)-1]*Size(3)/(size(Map,3)-1); %z node positions

%Plot side view of route
subplot(3,3,3)
plot([0;gd],Path(:,3),'r.-','LineWidth',1.5,'MarkerSize',13)
xlabel('Ground Distance (km)','fontsize',14)
ylabel('z (km)','fontsize',14)
grid minor
hold off

%Plot plan view of route
subplot(3,3,[6 9])
plot(Path(:,1),Path(:,2),'r.-','LineWidth',1.5,'MarkerSize',13)
hold on
contourslice(x,y,z,Map,[],[],0) 
contourslice(x,y,z,Map,[],[],0.2) 
contourslice(x,y,z,Map,[],[],0.4) 
contourslice(x,y,z,Map,[],[],0.6) 
contourslice(x,y,z,Map,[],[],0.8) 
contourslice(x,y,z,Map,[],[],1) 
xlabel('x (km)','fontsize',14)
ylabel('y (km)','fontsize',14)
axis([0 Size(1) 0 Size(2)])
daspect([20 20 1])
grid minor
hold off
set(gcf,'position',[10,10,1000,500])
end