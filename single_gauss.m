s = 0;
mean = [double(0); double(0); double(0)];
files = ["106"; "114"; "121"; "137"; "144"; "152"; "160"; "168"; "176"; "192"; "200"; "208"; "216"; "223"; "231"; "248"; "256"; "264"; "280"; "68"; "76"; "91"];
test_file = "train_images\231.jpg";
% 68 -> 1823
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
                area = area + 1;
                mean = mean + [double(R(i,j)) ; double(G(i,j)) ; double(B(i,j))];
            end
        end
    end
    s = s + area;
    dists(filenum) = dist;
    areas(filenum) = area;
end

mu = mean / s;
cova = [0 0 0; 0 0 0; 0 0 0];

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
                x = [double(R(i,j)) ; double(G(i,j)) ; double(B(i,j))];
                cova = cova + (((x - mu) * transpose(x - mu)) / s);
            end
        end
    end
end

det_cova = det(cova);

rgbs = imread(test_file);
R = rgbs(:,:,1); 
G = rgbs(:,:,2); 
B = rgbs(:,:,3);
count = 0;
max_like = 0;
under = 0;
heatmap = []
plz_god = 0;
for i = 1:size(R)
    for j = 1:size(R,2)
        x = [double(R(i, j)) ; double(G(i, j)) ; double(B(i, j))];
        like = (10^6) * 2 * (1 / sqrt((2 * pi)^3 * det_cova)) * exp(-1 * 0.5 * transpose(x - mu) * inv(cova) * (x - mu));
        if like > max_like
            max_like = like
        end
        if like > 1
            count = count + 1;
        end
        heatmap(i,j) = like;
        plz_god = plz_god + 1;
    end
end
under
max_like
count
%imagesc(heatmap)
scatter(areas, dists)

