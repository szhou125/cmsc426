clear all
load('GMM.mat')
test_file = 'test_images/8.jpg';
rgbs = imread(test_file);
R = rgbs(:,:,1); 
G = rgbs(:,:,2); 
B = rgbs(:,:,3);
thresh = 0.7;
heatmap = [];
for i = 1:size(R)
    for j = 1:size(R,2)
        x = [double(R(i, j)) ; double(G(i, j)) ; double(B(i, j))];
        like = 0;
        for tmp = 1:k
            local_mu = [mus(tmp,1) ; mus(tmp,2) ; mus(tmp,3)];
            local_cova = [covars(tmp,1) covars(tmp,2) covars(tmp,3) ; covars(tmp,2) covars(tmp,4) covars(tmp,5) ; covars(tmp,3) covars(tmp,5) covars(tmp,6)];
            local_pi = pis(tmp);
            like = like + local_pi * (10^6) * (1 / sqrt((2 * pi)^3 * det(local_cova))) * exp(-1 * 0.5 * transpose(x - local_mu) * ((local_cova) \ (x - local_mu)));
        end
        heatmap(i,j) = like;
    end
end
binary_heatmap = heatmap > thresh;
props_table = regionprops('struct',binary_heatmap,'Centroid','MajorAxisLength','MinorAxisLength');
max_area = 0;
for i = 1:size(props_table, 1)
    local_area = pi * (props_table(i).MajorAxisLength) * (props_table(i).MinorAxisLength) / 4;
    if local_area > max_area
        max_area = local_area;
    end
end
f=fit(areas', dists', 'linearinterp');
f(max_area)
imagesc(binary_heatmap)
%imagesc(heatmap)