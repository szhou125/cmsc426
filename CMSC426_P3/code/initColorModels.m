function ColorModels = initializeColorModels(IMG, Mask, MaskOutline, LocalWindows, BoundaryWidth, WindowWidth)
% INITIALIZAECOLORMODELS Initialize color models.  ColorModels is a struct you should define yourself.
%
% Must define a field ColorowModels.Confidences: a cell array of the color confidence map for each local window.
    %is this even used
    lab = rgb2lab(IMG);
    MaskOutline_b = imdilate(MaskOutline, strel('disk',BoundaryWidth));

    
    numLocalWindows = size(LocalWindows, 1);
    confidence = cell(numLocalWindows, 1);
    b_model = cell(numLocalWindows, 1);
    f_model = cell(numLocalWindows, 1);
    b_points = cell(numLocalWindows, 1);
    f_points = cell(numLocalWindows, 1);
    distances = cell(numLocalWindows, 1);
    probs = cell(numLocalWindows, 1);
    
    for i = 1:numLocalWindows
  
        %gets corners of window
        lowerX = LocalWindows(i,1) - WindowWidth / 2;
        upperX = LocalWindows(i,1) + WindowWidth / 2;
        lowerY = LocalWindows(i,2) - WindowWidth / 2;
        upperY = LocalWindows(i,2) + WindowWidth / 2; 
        
        %account for boundary width?
        %get pixels for gmm
        f_pixels = [];
        b_pixels = [];
        for j = lowerX:upperX
            for k = lowerY:upperY
                if (Mask(k,j) == 1 && MaskOutline_b(k,j) == 0)
                    f_pixels = [f_pixels; lab(k,j,1), lab(k,j, 2), lab(k,j, 3)];
                elseif (Mask(k,j) == 0 && MaskOutline_b(k,j) == 0) 
                    b_pixels = [b_pixels; lab(k,j,1), lab(k,j, 2), lab(k,j, 3)];
                end
            end
        end
        
        %create gmm models for foreground and background
       
        f_gmm = fitgmdist(f_pixels, 1);
        b_gmm = fitgmdist(b_pixels, 1);

        %foreground probability matrix
       p_matrix = zeros(WindowWidth, WindowWidth);
       col_index = 1;
        for j = lowerX:upperX
            row_index = 1;
            for k = lowerY:upperY
                b_prob = pdf(b_gmm, [lab(k,j, 1), lab(k,j, 2), lab(k,j, 3)]);
                f_prob = pdf(f_gmm, [lab(k,j, 1), lab(k,j, 2), lab(k,j, 3)]);
                prob = f_prob / (f_prob + b_prob);
                p_matrix(row_index,col_index) = prob;
                row_index = row_index + 1;
            end
            col_index = col_index + 1;        
        end
        
        %distance to mask outline
        omega = bwdist(MaskOutline);
        omega = omega(lowerY:upperY, lowerX:upperX);
        distance = -omega.^2;
        omega = exp((distance) / ((WindowWidth / 2) ^ 2));
        
        
        %calculate confidence
        fc = MaskOutline(lowerY:upperY,lowerX:upperX) - p_matrix;
        fc = abs(fc).*omega;
        fc = sum(fc(:));
        fc = 1 - fc / sum(omega(:));
        
        
        %adds to struct
        confidence{i} = fc;
        f_model{i} = f_gmm;
        b_model{i} = b_gmm;
        b_points{i} = b_pixels;
        f_points{i} = f_pixels;
        distances{i} = distance;
        probs{i} = p_matrix;
    end
    
    ColorModels = struct;
    ColorModels.Confidences = confidence;
    ColorModels.bGMM = b_model;
    ColorModels.fGMM = f_model;
    ColorModels.Windows = LocalWindows;
    ColorModels.fPoints = f_points;
    ColorModels.bPoints = b_points;
    ColorModels.distance = distances;
    ColorModels.prob = probs
end

