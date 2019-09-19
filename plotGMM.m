clear all
load('GMM.mat')
for i = 1:size(mus)
    local_cova = [covars(i, 1) covars(i, 2) covars(i, 3) ; covars(i, 2) covars(i, 4) covars(i, 5) ;  covars(i, 3) covars(i, 5) covars(i, 6) ];
    local_mu = [mus(i, 1) ; mus(i, 2) ; mus(i, 3)];
    [x,y,z] = ellipsoid(local_mu(1),local_mu(2),local_mu(3),sqrt(covars(i,1)),sqrt(covars(i,4)),sqrt(covars(i,6)));
    surf(x, y, z);
    hold on
end
hold off