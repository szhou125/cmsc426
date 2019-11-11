function [x, y] = equidistantPointsOnPerimeter(x_vertices, y_vertices, num_points)
% EQUIDISTANTPOINTSONPERIMETER Find equally spaced points along the perimeter of a polygon.
%
% Source: https://blogs.mathworks.com/steve/2012/07/06/walking-along-a-path/#5d684ab6-0edb-4c52-8df2-4bce4a79ee80

xy = [x_vertices y_vertices];
d = diff(xy,1);
dist_from_vertex_to_vertex = hypot(d(:,1), d(:,2));
cumulative_dist_along_path = [0; cumsum(dist_from_vertex_to_vertex,1)];
dist_steps = linspace(0, cumulative_dist_along_path(end), num_points);
points = interp1(cumulative_dist_along_path, xy, dist_steps);

x = points(:,1);
y = points(:,2);

end

