%% Load some nice data
[imTL,im410_1, im410_2,movement] = loadErrorData();

%% Demonstration of Segmentation

% Segmentation using static threshold
t = 2000;
i = 2;
I = im410_1(:,:,i);
mask = I;
mask(mask < 2000) = 0;
mask(mask <= 2000) = 1;
imshow(mask, []);