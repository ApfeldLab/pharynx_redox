imgs = loadTiffImageStack("/Users/sean/code/wormAnalysis/data/NEW_DATA/2017-02-22_HD233_SAY47/170222 HD233 SAY47.tif");
metadata = jsondecode(fileread("/Users/sean/code/wormAnalysis/data/NEW_DATA/2017-02-22_HD233_SAY47/metadata.json"));

nFrames = size(imgs, 3);
nFramesPerAnimal = size(metadata.lambda_sequence, 1);

% Ensure that the number of frames is divisible by the number of frames 
% specified by the imaging scheme in the JSON metadata
assert(mod(size(imgs, 3), nFramesPerAnimal) == 0, ...
    sprintf("The number of frames in the image stack (%u) is not divisible the number of frames (%u) described by the imaging scheme (%s).", nFrames, nFramesPerAnimal, string(join(metadata.lambda_sequence, "/"))));

% Ensure that there are no duplicates in the lambda sequence
assert(size(unique(metadata.lambda_sequence), 1) == size(metadata.lambda_sequence, 1), ...
    "Ensure that there are no duplicates in the lambda sequence (duplicates must be ordered e.g. 470_1, 410_1, 470_2, 410_2)")

image_stacks = struct;
for i=1:size(metadata.lambda_sequence,1)
    wavelength = strcat('i', metadata.lambda_sequence{i});
    image_stacks.(wavelength) = imgs(:, :, i:nFramesPerAnimal:end);
end