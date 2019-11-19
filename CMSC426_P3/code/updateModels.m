function [mask, LocalWindows, ColorModels, ShapeConfidences] = ...
    updateModels(...
        NewLocalWindows, ...
        LocalWindows, ...
        CurrentFrame, ...
        warpedMask, ...
        warpedMaskOutline, ...
        WindowWidth, ...
        ColorModels, ...
        ShapeConfidences, ...
        ProbMaskThreshold, ...
        fcutoff, ...
        SigmaMin, ...
        R, ...
        A ...
    )
% UPDATEMODELS: update shape and color models, and apply the result to generate a new mask.
% Feel free to redefine this as several different functions if you prefer.

    %update color models
    %%%%%%%%%%%%%%%%%%%%FIX X and Y IN ALL FILES!!!!
    numLocalWindows = size(NewLocalWindows, 1);
    confidence = cell(numLocalWindows, 1);
    foreground = cell(numLocalWindows, 1);
    background = cell(numLocalWindows, 1);
    f_prob = cell(numLocalWindows, 1);
    b_prob = cell(numLocalWindows, 1);
    distances = cell(numLocalWindows, 1);
    probs = cell(numLocalWindows, 1);
    w = cell(numLocalWindows, 1);

    %need to update cells
    %can i condense this since its in 1 file
    %what are the fields of each?
    
    img = rgb2lab(CurrentFrame);

    for i = 1:size(NewLocalWindows, 1)
        
        
        %gets corners of window
        lowerX = round(LocalWindows(i,1) - WindowWidth / 2);
        upperX = round(LocalWindows(i,1) + WindowWidth / 2);
        lowerY = round(LocalWindows(i,2) - WindowWidth / 2);
        upperY = round(LocalWindows(i,2) + WindowWidth / 2); 
        
        window = img(lowerY:upperY,lowerX:upperX,:);
        
        l = window(:,:,1);
        a = window(:,:,2);
        b = window(:,:,3);

        restructuredPoints(:,1:3) = [l(:),a(:),b(:)];
        
        %for each point in restructured points
        curr_window_f_prob = pdf(ColorModels.fGMM{i}, restructuredPoints);
        curr_window_b_prob = pdf(ColorModels.bGMM{i}, restructuredPoints);
        
        
        old_probs = curr_window_f_prob ./ (curr_window_f_prob + curr_window_b_prob);
        
        f_indices = find(old_probs >= 0.75);
        b_indices = find(old_probs <= 0.25);
        
        %dimensions of probs?
        f_probs = restructuredPoints(f_indices,:);
        b_probs = restructuredPoints(b_indices,:);
        
        all_f = [f_probs; ColorModels.fPoints{i}];
        all_b = [b_probs; ColorModels.bPoints{i}];
        
        fGMM = fitgmdist(all_f, 1);        
        bGMM = fitgmdist(all_b, 1);
        
        new_window_f_prob = pdf(fGMM, restructuredPoints);
        new_window_b_prob = pdf(bGMM, restructuredPoints);
        new_probs = new_window_f_prob ./ (new_window_f_prob + new_window_b_prob);

        %just gets indices
        %should this be 1?
        %probability is 1?
        new_f_count = find(new_probs == 1);
        %can i rewrite this lol like total - f count or something?
        new_b_count = find(new_probs ~= 1);  
        
        NewDistance = bwdist(warpedMaskOutline);
        NewDistance = NewDistance( lowerY:upperY, lowerX:upperX);
        NewDistance = -NewDistance.^2;
        
        if(size(new_f_count,1) < 1.2 * size(ColorModels.fPoints{i}, 1))
            newfGMM = fGMM;
            newbGMM = bGMM;
            newFPoints = restructuredPoints(new_f_count, :);
            newBPoints = restructuredPoints(new_b_count, :);
            
            omega = exp(NewDistance / ((WindowWidth/2)^2));
            omega = omega(:);
            
            conf = warpedMaskOutline( lowerY:upperY, lowerX:upperX);
            conf = conf(:);
            conf = abs( conf - new_probs ).*omega;
            conf = (1 - sum( conf(:))) / (sum(omega(:)));
            
            %do i need to use size?
            prob = reshape(new_probs, [size(window, 1) size(window, 2)]);
            
        else
            %set confidence, fgmm, bgmm, fpx, bpx, newprob
            newfGMM = ColorModels.fGMM{i};
            newbGMM = ColorModels.bGMM{i};
            conf = ColorModels.Confidences{i};
            newFPoints = f_probs;
            newBPoints = b_probs
            prob = reshape(old_probs, [size(window, 1) size(window, 2)]);
        end
        w{i} = NewLocalWindows(i,:);
        confidence{i} = conf;
        foreground{i} = newFPoints;
        background{i} = newBPoints;
        f_prob{i} = newfGMM;
        b_prob{i} = newbGMM;
        distances{i} = NewDistance;
        probs{i} = prob;
    end
    
    newColorModel = struct;
    newColorModel.Confidences = confidence;
    newColorModel.bGMM = b_prob;
    newColorModel.fGMM = f_prob;
    newColorModel.Windows = w;
    newColorModel.fPoints = foreground;
    newColorModel.bPoints = background;
    newColorModel.distance = distances;
    newColorModel.prob = probs;

    %uses the new color models
    ShapeConfidences = initShapeConfidences(NewLocalWindows, newColorModel, WindowWidth, SigmaMin, A, R, fcutoff);
    
    %is the x and y switched?
    %get warped mask windows
    LocalWarpedMasks = cell(numLocalWindows, 1);
    for i = 1:numLocalWindows
        lowerX = round(NewLocalWindows(i,1) - WindowWidth / 2);
        upperX = round(NewLocalWindows(i,1) + WindowWidth / 2);
        lowerY = round(NewLocalWindows(i,2) - WindowWidth / 2);
        upperY = round(NewLocalWindows(i,2) + WindowWidth / 2); 
        LocalWarpedMasks{i} = warpedMask( lowerY:upperY, lowerX:upperX);
    end
    
    %merge models
    pfs = cell(numLocalWindows, 1);
    for i = 1:numLocalWindows
        localMask = LocalWarpedMasks{i};
        shapeConf = ShapeConfidences.Confidence{i}
        %what does colormodels have?
        pfs{i} = shapeConf.*localMask + (1 - shapeConf) .* newColorModel.prob{i};
    end
    
    %merging
    maskNum = zeros(size(warpedMask));
    maskDenom = zeros(size(warpedMask));
    for i = 1:numLocalWindows
        pfsk = pfs{i};
        
        lowerX = round(NewLocalWindows(i,1) - WindowWidth / 2);
        upperX = round(NewLocalWindows(i,1) + WindowWidth / 2);
        lowerY = round(NewLocalWindows(i,2) - WindowWidth / 2);
        upperY = round(NewLocalWindows(i,2) + WindowWidth / 2); 
        
        for row = lowerY:upperY
            for col = lowerX:upperX
                
                %for j = 1:numLocalWindows
                %    checklowerX = round(NewLocalWindows(i,1) - WindowWidth / 2);
                %    checkupperX = round(NewLocalWindows(i,1) + WindowWidth / 2);
                %    checklowerY = round(NewLocalWindows(i,2) - WindowWidth / 2);
                %    checkupperY = round(NewLocalWindows(i,2) + WindowWidth / 2); 
                    
                %    if (NewLocalWindows(i,1) <= checkupperX &&...
                %            NewLocalWindows(i,1) >= checklowerX &&...
                %            NewLocalWindows(i,2) <= checkupperY &&...
                %            NewLocalWindows(i,2) >= checklowerY) 
                        denom = inv(pdist([NewLocalWindows(i,1), NewLocalWindows(i,2); col, row]) + 0.1);
                        %switch x and y?
                        num = pfsk( col - lowerX + 1, row - lowerY + 1) * denom;
                        maskNum(row, col) = maskNum(row,col) + num;
                        maskDenom(row, col) = maskDenom(row,col) + denom;
                %    end
                        
                    
                %end
                

            end
        end  
    end
    
    maskDenom(maskDenom == 0) = 1;
    Pfs = maskNum ./ maskDenom;
    
    %copy mask
    %update lowerX/upperX/lowerY/upperY
    UpdatedMask = zeros(size(Pfs));
    UpdatedMask(Pfs > ProbMaskThreshold) = 1;
    for i = 1:size(LocalWindows, 1)
        lowerY = round(LocalWindows(i,2) - WindowWidth / 2);
        upperY = round(LocalWindows(i,2) + WindowWidth / 2);
        lowerX = round(LocalWindows(i,1) - WindowWidth / 2);
        upperX = round(LocalWindows(i,1) + WindowWidth / 2); 
        window = UpdatedMask(lowerY:upperY, lowerX:upperX);
        window = imfill(window, 'holes');
        UpdatedMask(lowerY:upperY, lowerX:upperX) = window;
        
    end
    mask = imfill(UpdatedMask, 'holes');  
    LocalWindows = NewLocalWindows;
    ColorModels = newColorModel;
    
end

