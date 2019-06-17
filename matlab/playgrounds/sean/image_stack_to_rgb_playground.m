outputVideo = VideoWriter('~/Desktop/imE.avi');
outputVideo.FrameRate = 3;
open(outputVideo);
for i=1:size(imE_rgb, 4)
    writeVideo(outputVideo, imE_rgb(:,:,:,i));
end
close(outputVideo);

%%
imE = load('../data/tmp/imE.mat');
imE = imE.imE;

%%
map = colorcet('D6');
imE_rgb = imresize(im2uint8(imageStackToRGB(imE, map)), 2);
outputFileName = '~/Desktop/imE.tif';
options.color = true;
saveastiff(imE_rgb, outputFileName, options);