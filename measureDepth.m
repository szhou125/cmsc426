area = 500; % change to whatever area (in pixels) that you want to find distance for
distances = [106 114 121 137 144 152 160 168 176 192 200 208 216 223 231 248 256 264 280 68 76 91 99];
total_pixels = [];

for i = 1:23
   load("data\" + distances(i) + ".mat");
   total = sum(b(:));
   total_pixels(i) = total;
   plot(total, distances(i), "ok");
   hold on
end

f=fit(transpose(total_pixels), transpose(distances), 'linearinterp');
plot(f);

f(area)

hold off