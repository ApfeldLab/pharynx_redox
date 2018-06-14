%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration, global variables
N_CHANNELS = 3;

CHANNEL_TL  = 1;
CHANNEL_470 = 2;
CHANNEL_410 = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing

images = tiffread2('/Users/sean/Desktop/matlab_2018_06_07_SAY98_HD233/2018_06_07_SAY98_HD233.tif');
images = cat(3,images.data);
% images is a matrix of [height,width,n], where n is the number of images
% so e.g. to display the 42nd image from a stack called I, try
% imshow(I(:,:,42);


% Separate into stacks
iTL  = images(:, :,  CHANNEL_TL:N_CHANNELS:end);
i410 = images(:, :, CHANNEL_410:N_CHANNELS:end);
i470 = images(:, :, CHANNEL_470:N_CHANNELS:end);

