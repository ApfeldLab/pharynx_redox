
% Get channels from user
channel_list = input('Enter the image channels in order, separated by commas (default is [TL, 410, 470]): ', 's');
if strlength(channel_list) == 0
    channel_list = 'TL,410,470';
end
channel_list = strsplit(channel_list, ',');
% TODO: save this list to directory

% Load raw image stack
directory = '/Users/sean/Desktop/vab1_2do_05_25_mvmt/';
image_name = 'vab1_2do_05_25_mvmt.tif';
raw_stack = tiffread2(strcat(directory, image_name));
n_images = size(raw_stack, 2);
n_channels = size(channel_list, 2);

if mod(n_images, n_channels) ~= 0
    Error('Number of images in stack is not divisible by number of channels');
end