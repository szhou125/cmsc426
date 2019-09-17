s = 0;
mean = [double(0); double(0); double(0)];
no_s = 0;
no_mean = [double(0); double(0); double(0)];
prior = 0.99999;
files = ["106"; "114"; "121"; "137"; "144"; "152"; "160"; "168"; "176"; "192"; "200"; "208"; "216"; "223"; "231"; "248"; "256"; "264"; "280"; "68"; "76"; "91"];
test_file = "train_images\68.jpg";
dists = [];
areas = [];
for filenum = 1:size(files)
    filename = append('data\', files(filenum), '.mat');
    load(filename);
    imgname = 'train_images\' + files(filenum) + '.jpg';
    rgbs = imread(imgname);
    dist = str2num(files(filenum));
    area = 0;
    R = rgbs(:,:,1); 
    G = rgbs(:,:,2); 
    B = rgbs(:,:,3);
    
    for i = 1:size(b,1)
        for j = 1:size(b,2)
            if b(i,j) > 0
                s = s + 1;
                area = area + 1;
                mean = mean + [double(R(i,j)) ; double(G(i,j)) ; double(B(i,j))];
            else
                no_s = no_s + 1;
                no_mean = no_mean + [double(R(i,j)) ; double(G(i,j)) ; double(B(i,j))];
            end
        end
    end
    dists(filenum) = dist;
    areas(filenum) = area;
end
mu = [double(mean(1)); double(mean(2)); double(mean(3))] / s;
no_mu = [double(no_mean(1)) ; double(no_mean(2)) ; double(no_mean(3))] / no_s;

diff_R_R = double(0);
diff_R_G = double(0);
diff_R_B = double(0);
diff_G_R = double(0);
diff_G_G = double(0);
diff_G_B = double(0);
diff_B_R = double(0);
diff_B_G = double(0);
diff_B_B = double(0);

no_diff_R_R = double(0);
no_diff_R_G = double(0);
no_diff_R_B = double(0);
no_diff_G_R = double(0);
no_diff_G_G = double(0);
no_diff_G_B = double(0);
no_diff_B_R = double(0);
no_diff_B_G = double(0);
no_diff_B_B = double(0);

for filenum = 1:size(files)
    filename = append('data\', files(filenum), '.mat');
    load(filename);
    imgname = 'train_images\' + files(filenum) + '.jpg';
    rgbs = imread(imgname);
    R = rgbs(:,:,1); 
    G = rgbs(:,:,2); 
    B = rgbs(:,:,3);
    
    for i = 1:size(b,1)
        for j = 1:size(b,2)
            if b(i,j) > 0
                diff_R_R = diff_R_R + (double(R(i,j)) - mu(1)) * (double(R(i,j)) - mu(1));
                diff_R_G = diff_R_G + (double(R(i,j)) - mu(1)) * (double(G(i,j)) - mu(2));
                diff_R_B = diff_R_B + (double(R(i,j)) - mu(1)) * (double(B(i,j)) - mu(3));
                diff_G_R = diff_G_R + (double(G(i,j)) - mu(2)) * (double(R(i,j)) - mu(1));
                diff_G_G = diff_G_G + (double(G(i,j)) - mu(2)) * (double(G(i,j)) - mu(2));
                diff_G_B = diff_G_B + (double(G(i,j)) - mu(2)) * (double(B(i,j)) - mu(3));
                diff_B_R = diff_B_R + (double(B(i,j)) - mu(3)) * (double(R(i,j)) - mu(1));
                diff_B_G = diff_B_G + (double(B(i,j)) - mu(3)) * (double(G(i,j)) - mu(2));
                diff_B_B = diff_B_B + (double(B(i,j)) - mu(3)) * (double(B(i,j)) - mu(3));
            else
                no_diff_R_R = no_diff_R_R + (double(R(i,j)) - mu(1)) * (double(R(i,j)) - mu(1));
                no_diff_R_G = no_diff_R_G + (double(R(i,j)) - mu(1)) * (double(G(i,j)) - mu(2));
                no_diff_R_B = no_diff_R_B + (double(R(i,j)) - mu(1)) * (double(B(i,j)) - mu(3));
                no_diff_G_R = no_diff_G_R + (double(G(i,j)) - mu(2)) * (double(R(i,j)) - mu(1));
                no_diff_G_G = no_diff_G_G + (double(G(i,j)) - mu(2)) * (double(G(i,j)) - mu(2));
                no_diff_G_B = no_diff_G_B + (double(G(i,j)) - mu(2)) * (double(B(i,j)) - mu(3));
                no_diff_B_R = no_diff_B_R + (double(B(i,j)) - mu(3)) * (double(R(i,j)) - mu(1));
                no_diff_B_G = no_diff_B_G + (double(B(i,j)) - mu(3)) * (double(G(i,j)) - mu(2));
                no_diff_B_B = no_diff_B_B + (double(B(i,j)) - mu(3)) * (double(B(i,j)) - mu(3));
            end
        end
    end
end

cova = [diff_R_R diff_R_G diff_R_B ; diff_G_R diff_G_G diff_G_B ; diff_B_R diff_B_G diff_B_B]
no_cova = [no_diff_R_R no_diff_R_G no_diff_R_B ; no_diff_G_R no_diff_G_G no_diff_G_B ; no_diff_B_R no_diff_B_G no_diff_B_B]

cova = cova / s
no_cova = no_cova / no_s
det_cova = det(cova);
det_no_cova = det(no_cova);

rgbs = imread(test_file);
R = rgbs(:,:,1); 
G = rgbs(:,:,2); 
B = rgbs(:,:,3);
count = 0;
other_count = 0;
max_like = 0;
max_nolike = 0;

for i = 1:size(R)
    x = [double(R(i)) ; double(G(i)) ; double(B(i))];
    like = (1 / sqrt((2 * pi)^3 * det_cova)) * exp(-1 * 0.5 * transpose(x - mu) * inv(cova) * (x - mu));
    no_like = (1 / sqrt((2 * pi)^3 * det_no_cova)) * exp(-1 * 0.5 * transpose(x - no_mu) * inv(no_cova) * (x - no_mu));
    if (like * prior) > (no_like * (1 - prior))
        count = count + 1;
    end
    if (like ) > (no_like)
        other_count = 1000;
    end
    if like > max_like
        max_like = like;
    end
    if no_like > max_nolike
        max_nolike = no_like;
    end
end
max_like
max_nolike

count

scatter(areas, dists)

