%code snippets that I used

%to see how accurate my mask was
turtle = bsxfun(@times, img1, cast(mask1,class(img1)));

%to visualize local windows with mask on image
%[maskoutline1, localsamples1] = initLocalWindows(img1,mask1,35, 30, true)

%to read and create mask for image
%img1 = imread("1.jpg");
%mask1 = roipoly(img1);

m = initializeColorModels(img1, mask1, maskoutline1, localsamples1, 5, 30)
