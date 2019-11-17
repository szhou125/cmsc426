 function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.

% TODO
    WarpedPrevFrame = rgb2gray(WarpedPrevFrame);
    CurrentFrame = rgb2gray(CurrentFrame);
    
    ofEstimate = opticalFlowHS;
    %what is this line doing?
    estimateFlow(ofEstimate, WarpedPrevFrame);
    flow = estimateFlow(ofEstimate, CurrentFrame);
    
    velocityX = flow.Vx;
    velocityY = flow.Vy;
    distances = bwdist(Mask);
    NewLocalWindows = [];
    for i = 1:size(LocalWindows,1)
        
        %gets corners of window
        lowerX = round(LocalWindows(i,1) - Width / 2);
        upperX = round(LocalWindows(i,1) + Width / 2);
        lowerY = round(LocalWindows(i,2) - Width / 2);
        upperY = round(LocalWindows(i,2) + Width / 2); 
        
        changeX = 0;
        changeY = 0;
        
        for j = lowerX:upperX
            for k = lowerY:upperY
                if (Mask(k,j) == 1)
                    changeX = changeX + velocityX(k,j);
                    changeY = changeY + velocityY(k,j);
                end
            end
        end
        
        averageChangeX = changeX / (Width * Width);
        averageChangeY = changeY / (Width * Width);
        
        %is it 2 1 or 1 2
        NewLocalWindows(i, 1:2) =  [round(LocalWindows(i,2) + averageChangeX),
                                    round(LocalWindows(i,1) + averageChangeY)];
        
    end

end

