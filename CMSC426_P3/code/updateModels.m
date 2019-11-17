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
    %FIX X and Y IN ALL FILES!!!!
    numLocalWindows = size(NewLocalWindows, 1);
    confidence = cell(numLocalWindows, 1);
    foreground = cell(numLocalWindows, 1);
    background = cell(numLocalWindows, 1);
    windows_prob = cell(numLocalWindows, 1);
    f_prob = cell(numLocalWindows, 1);
    b_prob = cell(numLocalWindows, 1);
    distance = cell(numLocalWindows, 1);
    prob = cell(numLocalWindows, 1);

    %need to update cells
    %can i condense this since its in 1 file
    %need to really understand what im doing lol
    %what are the fields of each?
    for i = 1:size(NewLocalWindows, 1)
        
        img = rgb2lab(CurrentFrame);
        
        %gets corners of window
        lowerX = round(LocalWindows(i,1) - WindowWidth / 2);
        upperX = round(LocalWindows(i,1) + WindowWidth / 2);
        lowerY = round(LocalWindows(i,2) - WindowWidth / 2);
        upperY = round(LocalWindows(i,2) + WindowWidth / 2); 
        
        window = img(lowerX:upperX, lowerY:upperY,:);
        
        l = window(:,:,1);
        a = window(:,:,2);
        b = window(:,:,3);

        restructuredPoints(:,1:3) = [l(:),a(:),b(:)];
        
        curr_window_f_prob = pdf(ColorModels.fGMM{i}, restructuredPoints);
        curr_window_b_prob = pdf(ColorModels.bGMM{i}, restructuredPoints);
        
        old_probs = curr_window_f_prob ./ (curr_window_f_prob + curr_window_b_prob);
        
        f_indices = find(old_probs >= 0.75);
        b_indices = find(old_probs <= 0.25);
        
        %dimensions of probs?
        f_probs = old_probs(f_indices,:);
        b_probs = old_probs(b_indices,:);
        
        all_f = [f_probs; ColorModels.f_probs{i}];
        all_b = [b_probs; ColorModels.b_probs{i}];
        
        fGMM = fitgmdist(all_f, 1);        
        bGMM = fitgmdist(all_b, 1);
        
        new_window_f_prob = pdf(fGMM, restructuredPoints);
        new_window_b_prob = pdf(bGMM, restructuredPoints);
        new_probs = new_window_f_prob ./ (new_window_f_prob + new_window_b_prob);

        %just gets indices
        new_f_count = find(old_probs == 1);
        %can i rewrite this lol like total - f count or something?
        new_b_count = find(old_probs ~= 1);  
        
        NewDistance = bwdist(MaskOutline);
        NewDistance = NewDistance(lowerX:upperX, lowerY:upperY);
        NewDistance = -NewDistance.^2;
        
        if(size(new_f_count,1) < 1.2 * size(ColorModels.f_probs{i}, 1))
            newfGMM = fGMM;
            newbGMM = bGMM;
            newFPoints = restructuredPoints(new_f_count, :);
            newBPoints = restructuredPoints(new_b_count, :);
            
            omega = exp(NewDistance / ((WindowWidth/2)^2));
            omega = omega(:);
            
            conf = MaskOutline(lowerX:upperX, lowerY:upperY);
            conf = conf(:);
            conf = abs( conf - new_probs ).*omega;
            conf = (1 - sum( conf(:))) / (sum(omega(:)));
            
            %do i need to use size?
            prob = reshape(new_probs, [size(window, 1) size(window, 2)]);
            
        else
            %set new confidence, fgmm, bgmm, fpx, bpx, newprob
            
        end
        
    end
    
    %uses the new color models
    ShapeConfidences = initShapeConfidences(NewLocalWindows, ColorModels, WindowWidth, SigmaMin, A, R, fcutoff);
    
    %is the x and y switched?
    LocalWarpedMasks = cell(numLocalWindows, 1);
    for i = 1:numLocalWindows
        lowerX = round(NewLocalWindows(i,1) - WindowWidth / 2);
        upperX = round(NewLocalWindows(i,1) + WindowWidth / 2);
        lowerY = round(NewLocalWindows(i,2) - WindowWidth / 2);
        upperY = round(NewLocalWindows(i,2) + WindowWidth / 2); 
        LocalWarpedMasks{i} = warpedMask(lowerX:upperX, lowerY:upperY);
    end
    
    pfs = cell(numLocalWindows, 1);
    for i = 1:numLocalWindows
        localMask = LocalWarpedMasks{i};
        shapeConf = ShapeConfidences.Confidence{i}
        %what does colormodels have?
        pfs{i} = shapeConf.*localMask + (1 - shapeConf) .* ColorModels.prob{i};
    end
    
    %merging
    for i = 1:numLocalWindows
        pfsk = pfs{i};
        
        lowerX = round(NewLocalWindows(i,1) - WindowWidth / 2);
        upperX = round(NewLocalWindows(i,1) + WindowWidth / 2);
        lowerY = round(NewLocalWindows(i,2) - WindowWidth / 2);
        upperY = round(NewLocalWindows(i,2) + WindowWidth / 2); 
        
        for row = lowerX:upperX
            for col = lowerY:upperY
                denom = inv(pdist([newLocalWindows(i,1), newLocalWindows(i,1); col, row]) + 0.1);
                %switch x and y?
                num = pkfs( j - lowerY, i - lowerX) * denom;
                maskNum(row, col) = maskNum(row, col) + num;
                maskDenom(row, col) = maskDenom(row, col) + denom;
            end
        end  
    end
    
    maskDenom(maskDenom == 0) = 1;
    Pfs = maskNum ./ maskDenom;
    
    %copy mask
    %update lowerX/upperX/lowerY/upperY
    UpdatedMask(Pfs > ProbMaskThreshold) = 1;
    for i = 1:size(LocalWindows, 1)
        lowerX = round(LocalWindows(i,2) - WindowWidth / 2);
        upperX = round(LocalWindows(i,2) + WindowWidth / 2);
        lowerY = round(LocalWindows(i,1) - WindowWidth / 2);
        upperY = round(LocalWindows(i,1) + WindowWidth / 2); 
        
        window = UpdatedMask(lowerX:upperX, lowerY:upperY);
        window = imfill(window, 'holes');
        UpdatedMask(lowerX:upperX, lowerY:upperY) = window;
        
    end
    mask = imfill(UpdatedMask, 'holes');    
    
end

