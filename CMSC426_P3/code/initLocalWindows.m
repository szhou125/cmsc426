function [MaskOutline, LocalSamples] = initLocalWindows(IMG,Mask,Num_windows,Window_width,ShowWindows)
% INITIALIZElOCALWINDOWS Create local windows on boundary of mask
%   
% If ShowWindows is true, shows a plot of the windows on the img.

MaskOutline = bwperim(Mask,4);

Boundaries = bwboundaries(Mask);

PointsForBoundary = zeros(length(Boundaries));
AllCoords = 0;

for i=1:length(Boundaries)
    AllCoords = AllCoords + length(Boundaries{i});
end

for i=1:length(Boundaries)
    Coords = length(Boundaries{i});
    PointsForBoundary(i) = round((Coords/AllCoords)*Num_windows);
    AllCoords = AllCoords - Coords;
    Num_windows = Num_windows - PointsForBoundary(i);
end

AllWindowCentersX = [];
AllWindowCentersY = [];

for i=1:length(Boundaries)
    Boundary = Boundaries{i};
    [WindowCentersX, WindowCentersY] = ...
        equidistantPointsOnPerimeter(Boundary(:,2), Boundary(:,1), PointsForBoundary(i));

    AllWindowCentersX = [AllWindowCentersX; AllWindowCentersX];
    AllWindowCentersY = [AllWindowCentersY; AllWindowCentersY];
end

LocalSamples = [AllWindowCentersX AllWindowCentersY];
LocalSamples = round(LocalSamples);

if ShowWindows 
    imshow(IMG);
    hold on
end
    
for i = 1:length(LocalSamples)
    x = LocalSamples(i,1);
    y = LocalSamples(i,2);
    yRange = (y-(Window_width/2)):(y+(Window_width/2 - 1));
    xRange = (x-(Window_width/2)):(x+(Window_width/2 - 1));
    
    if ShowWindows 
        plot(x, y, '.');
        rectangle('Position', [(x - Window_width/2) (y - Window_width/2) Window_width Window_width]);
    end
end

if ShowWindows
    hold off
end

end
