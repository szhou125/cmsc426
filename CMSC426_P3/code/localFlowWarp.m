function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.

% TODO
    WarpedPrevFrame = rgb2gray(WarpedPrevFrame);
    CurrentFrame = rgb2gray(CurrentFrame);
    
    ofEstimate = opticalFlowHS;
    estimateFlow(ofEstimate, WarpedPrevFrame);
    flow = estimateFlow(ofEstimate, CurrentFrame);
    
    

end

