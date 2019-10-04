function [pano] = MyPanorama()
    % Must load images from ../Images/Input/
    % Must return the finished panorama.

    % Curr image is image 1
    % TODO change image read path to ../Images/Input/
    % (It's set to this path for testing purposes
    location = "..\Images\train_images\Set1\";
    curr_img = imread(location + "1.jpg");

    a=dir(location + '*.jpg');
    num_images=size(a,1);

    for I = 2:num_images
        % merge next_img into curr_img
        next_img = imread(location + I + ".jpg");
        % Corner detection -- convert to grayscale
        curr_img_grayscale = rgb2gray(curr_img);
        next_img_grayscale = rgb2gray(next_img);
        % Detect corners on both images
        curr_img_corners = cornermetric(curr_img_grayscale);
        next_img_corners = cornermetric(next_img_grayscale);
        % Get best regional corners out of these on both images
        curr_img_maxes = imregionalmax(curr_img_corners);
        next_img_maxes = imregionalmax(next_img_corners);

        % Run ANMS
        [curr_img_x_best, curr_img_y_best] = ANMS(curr_img, curr_img_maxes, curr_img_corners);
        [next_img_x_best, next_img_y_best] = ANMS(next_img, next_img_maxes, next_img_corners);

        % Display ANMS results -- uncomment this line to do so
        %displayANMS(next_img, next_img_x_best, next_img_y_best)
        % Get feature descriptors
        curr_img_descriptors = feature_descriptors(curr_img_grayscale, curr_img_x_best, curr_img_y_best);
        next_img_descriptors = feature_descriptors(next_img_grayscale, next_img_x_best, next_img_y_best);
        % Now match those features with minimum SSD
        [matchedPoints1, matchedPoints2, best_matches] = feature_matching(curr_img_descriptors, next_img_descriptors, curr_img_x_best, curr_img_y_best, next_img_x_best, next_img_y_best);
        % Uncomment this line to see the matched features without RANSAC
        %showMatchedFeatures(curr_img, next_img, matchedPoints1, matchedPoints2, 'montage');

        % RANSAC to make it better
        [r_matchedPoints1, r_matchedPoints2, H_best] = RANSAC(matchedPoints1, matchedPoints2, curr_img_descriptors, next_img_descriptors, best_matches, curr_img_x_best, curr_img_y_best, next_img_x_best, next_img_y_best);
        % Uncomment the following line to see the matches with RANSAC
        showMatchedFeatures(curr_img, next_img, r_matchedPoints1, r_matchedPoints2, 'montage');
        
        % H_best is the homography between the two images
        % r_matchedPoints1 and r_matchedPoints2 are each nx2 vectors, where
        %   n is the number of matched points (x, y). r_matchedPoints1 has
        %   points from the "curr_img", and r_matchedPoints2 has points
        %   from the "next_img".
        
        % Img1 (curr_img) is the destination. Img2 (next_img) is the src
        % Now, we warp and blend
        %blendedImage = getBlendedImage(...)
        % Set curr_img to blendedImage and continue the loop
    end
end
%{
function blendedImage = getBlendedImage(curr_img, next_img, r_matchedPoints1, r_matchedPoints2, H_best)
    % Apply homography to the corners of the next_img
    [next_y_max, next_x_max, ~] = size(next_img);
    [next_x_corners, next_y_corners] = apply_homography(H_best, [0;0;next_x_max;next_x_max], [0;next_y_max;0;next_y_max]);
    % TODO...
end
%}
function [r_matchedPoints1, r_matchedPoints2, H_best] = RANSAC(matchedPoints1, matchedPoints2, img1_descriptors, img2_descriptors, best_matches, img1_x_best, img1_y_best, img2_x_best, img2_y_best)
    % best_matches object from feature_matching. Gives us the indices of
    % the respective descriptors
    [~, num_matches] = size(best_matches);
    max_iterations = 25000;  % Seems to converge after a point, but this is fast
    curr_iteration = 0;
    best_inliers = -1;
    threshold = 400;  % Tried many different values, this works well!

    % Get in proper format for apply_homography
    all_src_x = reshape(img2_x_best(best_matches(2, :)), [num_matches, 1]);
    all_src_y = reshape(img2_y_best(best_matches(2, :)), [num_matches, 1]);
    all_dest_x = reshape(img1_x_best(best_matches(1, :)), [num_matches, 1]);
    all_dest_y = reshape(img1_y_best(best_matches(1, :)), [num_matches, 1]);
    
    best_inliers_src_x = [];
    best_inliers_src_y = [];
    best_inliers_dest_x = [];
    best_inliers_dest_y = [];

    while curr_iteration < max_iterations
        curr_iteration = curr_iteration + 1;
        cols = randperm(size(best_matches,2),4);
        columns = best_matches(:, cols);
        % Image 1 is the destination
        % Image 2 is the source
        dest_x = img1_x_best(columns(1, :));
        dest_y = img1_y_best(columns(1, :));
        src_x = img2_x_best(columns(2, :));
        src_y = img2_y_best(columns(2, :));
        % Get homography
        H = est_homography(dest_x, dest_y, src_x, src_y);
        % Now, apply homography and get SSD
        [new_x, new_y] = apply_homography(H, all_src_x, all_src_y);
        [count, inliers_src_x, inliers_src_y, inliers_dest_x, inliers_dest_y] = get_ssd_points(new_x, new_y, all_dest_x, all_dest_y, all_src_x, all_src_y, threshold);
        if or(best_inliers == -1, count > best_inliers)
            best_inliers = count;
            best_inliers_src_x = reshape(inliers_src_x, [count, 1]);
            best_inliers_src_y = reshape(inliers_src_y, [count, 1]);
            best_inliers_dest_x = reshape(inliers_dest_x, [count, 1]);
            best_inliers_dest_y = reshape(inliers_dest_y, [count, 1]);
        end
    end

    [~, num_best_matches] = size(best_matches);
    
    % Re-calculate homography only with the inliers
    H_best = est_homography(best_inliers_dest_x, best_inliers_dest_y, best_inliers_src_x, best_inliers_src_y);
    % Now, re-apply it to those points
    [new_dest_x, new_dest_y] = apply_homography(H_best, best_inliers_src_x, best_inliers_src_y);
    
    r_matchedPoints2 = [best_inliers_src_x best_inliers_src_y];
    r_matchedPoints1 = [new_dest_x new_dest_y];
end

function [count, inliers_src_x, inliers_src_y, inliers_dest_x, inliers_dest_y] = get_ssd_points(src_x, src_y, all_dest_x, all_dest_y, all_src_x, all_src_y, thresh)
    [len, ~] = size(src_x);
    count = 0;
    inliers_src_x = [];
    inliers_src_y = [];
    inliers_dest_x = [];
    inliers_dest_y = [];
    for indx = 1:len
        x_diff = all_dest_x(indx) - src_x(indx);
        y_diff = all_dest_y(indx) - src_y(indx);
        dist = (x_diff * x_diff) + (y_diff * y_diff);
        if dist < thresh
            count = count + 1;
            inliers_src_x = [inliers_src_x all_src_x(indx)];
            inliers_src_y = [inliers_src_y all_src_y(indx)];
            inliers_dest_x = [inliers_dest_x all_dest_x(indx)];
            inliers_dest_y = [inliers_dest_y all_dest_y(indx)];
        end
    end
end

function [matchedPoints1, matchedPoints2, best_matches] = feature_matching(img1_descriptors, img2_descriptors, img1_x_best, img1_y_best, img2_x_best, img2_y_best)
    %X = I1 - I2;
    %ssd = sum(X(:).^2); -- to find SSD
    best_matches = [];
    % sortrows(best_matches.',3).' where first row is i, second is j, third is ssd
    [~, n_best] = size(img1_x_best);
    for i = 1:n_best  % i is the indx of the img1 descriptors
        best_j = -1;
        best_ssd = -1;
        second_best_j = -1;
        second_best_ssd = -1;
        if img1_descriptors(64, i) == 0 && img1_descriptors(37, i) == 0 && img1_descriptors(1, i) == 0
            % If descriptor is all zeros, skip it
            continue
        end
        desc1 = img1_descriptors(:, i);
        for j = 1:n_best  % j is the indx of the img2 descriptors
            if img2_descriptors(64, j) == 0 && img2_descriptors(37, j) == 0 && img2_descriptors(1, j) == 0
                % If descriptor is all zeros, skip it
                continue
            end
            desc2 = img2_descriptors(:, j);
            diff = double(desc1) - double(desc2);
            ssd = sum(diff(:).^2);
            if or(best_j == -1, ssd < best_ssd)
                second_best_ssd = best_ssd;
                second_best_j = best_j;
                best_ssd = ssd;
                best_j = j;
            elseif or(second_best_j == -1, ssd < second_best_ssd)
                second_best_ssd = ssd;
                second_best_j = j;
            end
        end
        
        % We want the distances to be far apart
        % We have our best j for this one, compare to second best
        if second_best_ssd - best_ssd > 0.5
            v = reshape([i best_j (second_best_ssd - best_ssd)], [3 1]);
            best_matches = [best_matches v];
        end
    end
    
    [~, num_best_matches] = size(best_matches);
    matchedPoints1 = zeros(num_best_matches, 2);
    matchedPoints2 = zeros(num_best_matches, 2);
    best_matches = sortrows(best_matches.',3).';
    % Format
    for indx = 1:num_best_matches
        match = best_matches(:, indx);
        img1_x = img1_x_best(match(1));
        img1_y = img1_y_best(match(1));
        img2_x = img2_x_best(match(2));
        img2_y = img2_y_best(match(2));
        % Now, add to matchedPoints vectors
        matchedPoints1(indx, 1) = img1_x;
        matchedPoints1(indx, 2) = img1_y;
        matchedPoints2(indx, 1) = img2_x;
        matchedPoints2(indx, 2) = img2_y;
    end
end

function [descriptors] = feature_descriptors(img_grayscale, x_best, y_best)
    [~, n_best] = size(x_best);
    [y_size, x_size] = size(img_grayscale);
    descriptors = [];
    for indx = 1:n_best
        % Patch is of size 41x41, so point is the actual center
        x = x_best(indx);
        y = y_best(indx);
        % TODO handle patches near the edge?
        if or(or(or(x <= 20, y<=20), x + 20 > x_size), y + 20 > y_size)
            v = zeros(64,1);
            descriptors = [descriptors v];
            continue;
        end
        patch = img_grayscale(y-20:y+20, x-20:x+20);
        % Now, we take the gaussian blur!
        sig = 10;
        blurred = imgaussfilt(patch, sig);
        % Resize
        resized = imresize(blurred, [8 8]);
        % Reshape
        reshaped = double(reshape(resized, [64, 1]));
        % Now we need to standardize
        std_dev = std(reshaped);
        mean_reshaped = mean(reshaped);
        standardized = reshaped - mean_reshaped;
        standardized = standardized / std_dev;
        descriptors = [descriptors standardized];
    end
end

function [x_best, y_best] = ANMS(img, regional_maxes, corner_scores)
    n_best = 250;  % fixed number
    n_strong = 0;
    [y_size, x_size, ~] = size(img);
    x_max = [];
    y_max = [];
    for x = 1:x_size
        for y = 1:y_size
            if regional_maxes(y,x) == true
                n_strong = n_strong + 1;
                x_max = [x_max x];
                y_max = [y_max y];  % array concatenation
            end
        end
    end
    
    r = Inf(1, n_strong);
    for i = 1:n_strong
        for j = 1:n_strong
            if corner_scores(y_max(j), x_max(j)) > corner_scores(y_max(i), x_max(i))
                ED = ((x_max(j) - x_max(i))^2) + ((y_max(j) - y_max(i))^2);
                if ED < r(i)
                    r(i) = ED;
                end
            end
        end
    end
    
    x_best = zeros(1, n_best);
    y_best = zeros(1, n_best);
    
    % Now I need to sort r in descending order and pick lowest n_best points
    [~, I] = sort(r, 'descend');
    for i = 1:n_best
        indx = I(i);
        x_best(i) = x_max(indx);
        y_best(i) = y_max(indx);
    end
end


function [] = displayANMS(img, x_best, y_best)
    imshow(img);
    hold on
    
    plot(x_best, y_best, 'r.');
    hold off
end