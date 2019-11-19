function ShapeConfidences = initShapeConfidences(LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.
    
    numLocalWindows = size(LocalWindows, 1);
    f_s = cell(numLocalWindows);

    for i = 1:numLocalWindows
        if (ColorConfidences.Confidences{i} > fcutoff)
            sigmaS = SigmaMin + A*(ColorConfidences.Confidences{i} - fcutoff).^R;
        else
            sigmaS = SigmaMin;
        end
        
        f_s{i} = 1 - exp(-ColorConfidences.distance{i}.^2 / sigmaS^2);
    end
    
    ShapeConfidences = struct;
    ShapeConfidences.Confidence = f_s;
    
end
