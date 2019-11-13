function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,Windows)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
    img1 = rgb2gray(IMG1);
    img2 = rgb2gray(IMG2);
    
    points1 = detectHarrisFeatures(img1);
    points2 = detectHarrisFeatures(img2);
    
    [features1, valid_points1] = extractFeatures(img1, points1);
    [features2, valid_points2] = extractFeatures(img2, points2);
   
    indexPairs = matchFeatures(features1, features2);
    
    matched_points1 = valid_points1(indexPairs(:, 1), :);
    matched_points2 = valid_points2(indexPairs(:, 2), :);
    
    %poor feature matching
    %figure; showMatchedFeatures(img1,img2,matched_points1,matched_points2);
    
    %isolate foreground
    f_points1 = [];
    f_points2 = [];
    
    for i = 1:matched_points1.Count
        matched_point1 = matched_points1.Location(i,:);
        matched_point2 = matched_points2.Location(i,:);
        
        if (Mask(round(matched_point1(:,2)), round(matched_point1(:,1))) == 1) 
            f_points1 = [f_points1; matched_point1];
            f_points2 = [f_points2; matched_point2];
        end
        
    end
        
    %metric is 0?
    matched_points1 = cornerPoints(f_points1);
    matched_points2 = cornerPoints(f_points2);
    
    transform = estimateGeometricTransform(matched_points1, matched_points2, 'affine');
    
    %figure; showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
    
    WarpedMask = imwarp(Mask, transform);
    WarpedMaskOutline = bwperim(WarpedMask,4);
    WarpedFrame = imwarp(IMG1, transform);
    
    
    [x, y] = transformPointsForward(transform, Windows(:,1), Windows(:, 2));
    WarpedLocalWindows = [x, y];
       
    
end

