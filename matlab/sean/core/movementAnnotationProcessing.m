dataDir = '/Users/sean/code/wormAnalysis/data/NEW_DATA/';

subDirs = dir(dataDir);

for i=1:length(subDirs)
    if isfolder(subDirs(i).folder) && startsWith(subDirs(i).name, "20")
        experiment_dir = fullfile(subDirs(i).folder, subDirs(i).name);
        fprintf("Analyzing %s\n", experiment_dir);
        
        processForMovementAnnotation(experiment_dir);
    end
end

%%
processForMovementAnnotation("/Users/sean/code/wormAnalysis/data/NEW_DATA/2017_10_12-gld_1_RNAi_SAY88_HD233")

%
function processForMovementAnnotation(imageFileDirectory)
    imageFilePath = dir(fullfile(imageFileDirectory, '*.tif'));
    imageFilePath = fullfile(imageFileDirectory, imageFilePath(1).name);
    [~, baseImageName, ~] = fileparts(imageFilePath);
    
    allImages = loadTiffImageStack(imageFilePath);

%     iTL    = allImages(:,:,1:5:end);
    i470_1 = allImages(:,:,2:5:end);
    i410_1 = allImages(:,:,3:5:end);
    i470_2 = allImages(:,:,4:5:end);
    i410_2 = allImages(:,:,5:5:end);

    ratio_pair1 = double(i410_1) ./ double(i470_1);
    ratio_pair2 = double(i410_2) ./ double(i470_2);
    
%     seg_i470_1 = segmentPharynx(i470_1, 0, 2000);
    seg_i410_1 = segmentPharynx(i410_1, 0, 2000);
%     seg_i470_2 = segmentPharynx(i470_2, 0, 2000);
    seg_i410_2 = segmentPharynx(i410_2, 0, 2000);

    rot_i470_1 = rotatePharynx(i470_1, seg_i410_1);
    rot_i410_1 = rotatePharynx(i410_1, seg_i410_1);
    rot_i470_2 = rotatePharynx(i470_2, seg_i410_2);
    rot_i410_2 = rotatePharynx(i410_2, seg_i410_2);
    
    rot_ratio_pair1 = double(rot_i410_1) ./ double(rot_i470_1);
    rot_ratio_pair2 = double(rot_i410_2) ./ double(rot_i470_2);
    
    w = size(rot_ratio_pair1, 2);
    h = size(rot_ratio_pair1, 1);
    new_w = 70;
    new_h = 40;

    B = (h/2)-(new_h/2);
    T = (h/2)+(new_h/2);
    L = (w/2)-(new_w/2);
    R = (w/2)+(new_w/2);

    rot_ratio_pair1 = rot_ratio_pair1(B:T,L:R,:);
    rot_ratio_pair2 = rot_ratio_pair2(B:T,L:R,:);
    
    % Create mean-normalized ratio images
    
    % First, find the means inside each pharynx
    means_ratio_pair1 = calcMeanInsidePharynxes(seg_i410_1, ratio_pair1);
    means_ratio_pair2 = calcMeanInsidePharynxes(seg_i410_2, ratio_pair2);
    
    % Then, subtract the mean from each frame
    norm_rot_ratio_pair1 = (rot_ratio_pair1 - means_ratio_pair1) + 1;
    norm_rot_ratio_pair2 = (rot_ratio_pair2 - means_ratio_pair2) + 1;

    norm_ratio_pair1 = (ratio_pair1 - means_ratio_pair1) + 1;
    norm_ratio_pair2 = (ratio_pair2 - means_ratio_pair2) + 1;
    
    % Create composites
    comp_norm_rot_ratios = cat(1, norm_rot_ratio_pair1, norm_rot_ratio_pair2);
    comp_norm_ratios = cat(1, norm_ratio_pair1, norm_ratio_pair2);
    
    % Write all images to disk

    % split raw

    % ratio

    % rotated ratio

    % mean-normalized
    center = 1;
    buffer = 0.6;
    norm_img_range = [center - buffer, center + buffer];
    
    fprintf("Writing composite images to disk\n\n");
    
    writeImageStackToTiffStack(im2uint8(mat2gray(comp_norm_rot_ratios, norm_img_range)), ...
        fullfile(imageFileDirectory, strcat(baseImageName, '-comp_norm_rot_ratio.tif')));
    
    writeImageStackToTiffStack(im2uint8(mat2gray(comp_norm_ratios, norm_img_range)), ...
        fullfile(imageFileDirectory, strcat(baseImageName, '-comp_norm_ratio.tif')));
    
    writeImageStackToTiffStack(im2uint8(mat2gray(norm_ratio_pair1, norm_img_range)), ...
        fullfile(imageFileDirectory, strcat(baseImageName, '-norm_ratio_pair1.tif')));
    
    writeImageStackToTiffStack(im2uint8(mat2gray(norm_ratio_pair2, norm_img_range)), ...
        fullfile(imageFileDirectory, strcat(baseImageName, '-norm_ratio_pair2.tif')));
    
    writeImageStackToTiffStack(im2uint8(mat2gray(norm_rot_ratio_pair1, norm_img_range)), ...
        fullfile(imageFileDirectory, strcat(baseImageName, '-norm_rot_ratio_pair1.tif')));
    
    writeImageStackToTiffStack(im2uint8(mat2gray(norm_rot_ratio_pair2, norm_img_range)), ...
        fullfile(imageFileDirectory, strcat(baseImageName, '-norm_rot_ratio_pair2.tif')));
end

% ip = implay(comp_norm_unnorm_ratio_pair1);
% ip.Visual.ColorMap.UserRange    = 1;
% ip.Visual.ColorMap.UserRangeMin = .5;
% ip.Visual.ColorMap.UserRangeMax = 1.5;

% Helper Functions
function means = calcMeanInsidePharynxes(segStack, imStack)
    mask = imStack;
    mask(~segStack) = NaN;
    means = mean(mask, [1 2], 'omitnan');
%     means = squeeze(means(1,1,:));
end

function writeImageStackToTiffStack(imgStack, filepath)
    for K=1:length(imgStack(1, 1, :))
       imwrite(imgStack(:, :, K), filepath, 'WriteMode', 'append',  'Compression','none');
    end
end