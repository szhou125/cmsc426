function [] = showColorConfidences(IMG,MaskOutline,ColorConfidences,LocalWindows,Width);
% SHOWCOLORCONFIDENCES: plots the color confidence for each local window.
%
% ColorConfidencs should be a cell array of color confidences (one per window).

[W,H,~] = size(IMG);
ccNumer = zeros(W,H);
ccDenom = zeros(W,H);

const = 0.1;
[imH,imW] = size(MaskOutline);
for i = 1:length(LocalWindows)
    x = LocalWindows(i,1);
    y = LocalWindows(i,2);
    x = round(x); y = round(y);
    yMin = max((y-(Width/2)),1);
    yMax = min((y+(Width/2 - 1)),imH);
    xMin = max((x-(Width/2)),1); 
    xMax = min((x+(Width/2 - 1)),imW);
    yRange = yMin:yMax;
    xRange = xMin:xMax;
    
    WindowHeight = length(yRange);
    WindowWidth = length(xRange);

    if (WindowHeight < Width) || (WindowWidth < Width)
        continue;
    end
    
    
    center = zeros(Width, Width);
    center(Width/2+1,Width/2+1) = 1;
    distFromCenter = bwdist(center);

    LocalWeight = (distFromCenter + const).^(-1);
    LocalNum = LocalWeight .* repmat(ColorConfidences{i},[Width Width]);

    ccNumer(yRange, xRange) = ccNumer(yRange, xRange) + LocalNum;
    ccDenom(yRange, xRange) = ccDenom(yRange, xRange) + LocalWeight;
end

CC = ccNumer ./ ccDenom;
CC(isnan(CC)) = 0;

DilatedMaskOutline = imdilate(MaskOutline, strel('disk', 7, 4));

imagesc(CC .* DilatedMaskOutline);

end

