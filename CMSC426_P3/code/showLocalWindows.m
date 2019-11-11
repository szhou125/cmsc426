function [] = showLocalWindows(LocalWindows,Width,Type)

for i = 1:length(LocalWindows)
    x = LocalWindows(i,1);
    y = LocalWindows(i,2);
    
    plot(x, y, Type);
    rectangle('Position', [(x - Width/2) (y - Width/2) Width Width]);
end


end

