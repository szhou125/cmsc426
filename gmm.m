s = 0;
mu = [ 0 ; 0 ; 0];
files = ["106"; "114"; "121"; "137"; "144"; "152"; "160"; "168"; "176"; "192"; "200"; "208"; "216"; "223"; "231"; "248"; "256"; "264"; "280"; "68"; "76"; "91"];
test_file = "train_images\76.jpg";
% 68 -> 1823
dists = [];
areas = [];
R_balls = [];
G_balls = [];
B_balls = [];
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
                area = area + 1;
                mu = mu + [ double(R(i,j)) ; double(G(i,j)) ; double(B(i,j)) ];
                R_balls(size(R_balls,2) + 1) = double(R(i,j));
                G_balls(size(G_balls,2) + 1) = double(G(i,j));
                B_balls(size(B_balls,2) + 1) = double(B(i,j));
            end
        end
    end
    s = s + area;
    dists(filenum) = dist;
    areas(filenum) = area;
end

mu = mu / s;

k = 7;
tau = 0.5;
mus = [];
covars = [];
pis = [];

% Randomized Initialization

for i = 1:k
    mus(i, 1) = 255 * rand(1,1)%mu(1) + (80 * rand(1,1) - 40);
    mus(i, 2) = 255 * rand(1,1)%mu(2) + (80 * rand(1,1) - 40);
    mus(i, 3) = 255 * rand(1,1)%mu(3) + (80 * rand(1,1) - 40);
    covar_val = 1000*rand(1,1) + 1;
    covars(i, 1) = covar_val;
    covars(i, 2) = 0;
    covars(i, 3) = 0;
    covars(i, 4) = covar_val;
    covars(i, 5) = 0;
    covars(i, 6) = covar_val;
    pis(i) = rand(1,1) * (10^(-3));
end

sum_diff = 100;
old_mus = [];
iter = 0

while sum_diff > tau
    
    for i = 1:k
        old_mus(i, 1) = mus(i, 1);
        old_mus(i, 2) = mus(i, 2);
        old_mus(i, 3) = mus(i, 3);
    end
    
    alphas = [];
    
    for i = 1:k
        local_cova = [covars(i, 1) covars(i, 2) covars(i, 3) ; covars(i, 2) covars(i, 4) covars(i, 5) ; covars(i, 3) covars(i, 5) covars(i, 6)];
        local_mu = [mus(i, 1) ; mus(i, 2) ; mus(i, 3)];
        for j = 1:size(R_balls,2)
            x = [ R_balls(j) ; G_balls(j) ; B_balls(j) ];
            alphas(i,j) = pis(i) * (1 / sqrt((2 * pi)^3 * det(local_cova))) * exp(-1 * 0.5 * transpose(x - local_mu) * (local_cova \ (x - local_mu)));
        end
    end
    
    for j = 1:size(R_balls,2)
        local_sum = 0;
        for i = 1:k
            local_sum = local_sum + alphas(i, j);
        end
        for i = 1:k
            alphas(i, j) = alphas(i, j) / local_sum;
        end
    end
    
    for i = 1:k
        
        sum_alphas = 0;
        for j = 1:size(R_balls,2)
            sum_alphas = sum_alphas + alphas(i,j);
        end
        
        mus(i,1) = 0;
        mus(i,2) = 0;
        mus(i,3) = 0;
        for j = 1:size(R_balls,2)
            mus(i,1) = mus(i,1) + alphas(i,j) * R_balls(j);
            mus(i,2) = mus(i,2) + alphas(i,j) * G_balls(j);
            mus(i,3) = mus(i,3) + alphas(i,j) * B_balls(j);
        end
        for tmp = 1:3
            mus(i,tmp) = mus(i,tmp) / sum_alphas;
        end
        
        local_cova = [ 0  0 0 ; 0 0 0 ; 0 0 0];
        local_mu = [mus(i,1) ; mus(i,2) ; mus(i,3)];
        for j = 1:size(R_balls,2)
            local_x = [R_balls(j) ; G_balls(j) ; B_balls(j)];
            local_cova = local_cova + alphas(i,j) * (local_x - mu) * transpose(local_x - mu);
        end
        covars(i,1) = local_cova(1) / sum_alphas;
        covars(i,2) = local_cova(2) / sum_alphas;
        covars(i,3) = local_cova(3) / sum_alphas;
        covars(i,4) = local_cova(5) / sum_alphas;
        covars(i,5) = local_cova(6) / sum_alphas;
        covars(i,6) = local_cova(9) / sum_alphas;
        
        pis(i) = sum_alphas / size(R_balls,2);
        
    end
    
    sum_diff = 0;
    for i = 1:k
        new_full_mus = [mus(i,1) ; mus(i,2) ; mus(i,3)];
        old_full_mus = [old_mus(i,1) ; old_mus(i,2) ; old_mus(i,3)];
        sum_diff = sum_diff + norm(new_full_mus - old_full_mus);
    end
    
    iter = iter + 1
    
end

rgbs = imread("train_images\106.jpg");
R = rgbs(:,:,1); 
G = rgbs(:,:,2); 
B = rgbs(:,:,3);
count = 0;
thresh = 1;
heatmap = [];
for i = 1:size(R)
    for j = 1:size(R,2)
        x = [double(R(i, j)) ; double(G(i, j)) ; double(B(i, j))];
        like = 0;
        for tmp = 1:k
            local_mu = [mus(tmp,1) ; mus(tmp,2) ; mus(tmp,3)];
            local_cova = [covars(tmp,1) covars(tmp,2) covars(tmp,3) ; covars(tmp,2) covars(tmp,4) covars(tmp,5) ; covars(tmp,3) covars(tmp,5) covars(tmp,6)];
            local_pi = pis(tmp);
            like = like + local_pi * (10^8) * (1 / sqrt((2 * pi)^3 * det(local_cova))) * exp(-1 * 0.5 * transpose(x - local_mu) * inv(local_cova) * (x - local_mu));
        end
        if like > thresh
            count = count + 1;
        end
        heatmap(i,j) = like;
    end
end
%binary_heatmap = heatmap > 1;
%regionprops('table',binary_heatmap,'Centroid','MajorAxisLength','MinorAxisLength')
imagesc(heatmap)
%scatter(areas, dists)
