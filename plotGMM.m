clear all
load('GMM.mat')
for i = 1:size(mus)
    [x,y,z] = ellipsoid(mus(i, 1),mus(i, 2),mus(i, 3),sqrt(covars(i,1)),sqrt(covars(i,4)),sqrt(covars(i,6)));
    surf(x, y, z);
    hold on
end
xlabel('R')
ylabel('G')
zlabel('B')
hold off