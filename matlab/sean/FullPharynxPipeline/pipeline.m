%% Load Images
rawImgFilePath = "/Users/sean/Desktop/Analyses for Jodie's new Data/2018_10_14_SAY47_HD233_ctlRNAi_ifg1RNAi/2018_10_14_SAY47_HD233_ctlRNAi_ifg1RNAi.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);
imTable = struct2table(tiffread2(rawImgFilePath));
allImgs = cat(3, imTable.data{:}); % dimensions: WxHxi; e.g. first image is ims(:, :, 1)

imTL = allImgs(:, :, 1:3:end);
im410_raw = allImgs(:, :, 2:3:end);
im470_raw = allImgs(:, :, 3:3:end);
nAnimals = size(im410_raw, 3);

%% Subtract medians
median410 = median(im410_raw, [1 2]);
median470 = median(im470_raw, [1 2]);

im410 = im410_raw;
im470 = im470_raw;
% TODO vectorize
for i=1:nAnimals
    im410(:,:,i) = im410(:,:,i) - median410(i);
    im470(:,:,i) = im470(:,:,i) - median470(i);
end

%% Strain Names
strainTable = readtable(fullfile(rawImgDir, 'strains.csv'));
strains = cell(1,nAnimals);
for i=1:size(strainTable, 1)
    row_ = strainTable(i,:);
    strains(row_.Start:row_.End) = row_.Strain;
end
strains = categorical(cellstr(strains.'));
[GN, ~, G] = unique(strains);

%% Process Images

% This segmentation usually works, if the data looks a little wonky,
% take a look at the segmentation. That's usually the problem. If the
% segmentation fails. Go to ImageJ and try segmenting there, then load the
% masks here instead of calling segmentPharynx.
seg410 = segmentPharynx(im410, 0);
seg470 = segmentPharynx(im470, 0);

% [centroids, orientations] = calcCentroidsAndOrientations(seg410);
% 
% rotated410 = rotatePharynx(im410, centroids, orientations);
% rotated470 = rotatePharynx(im470, centroids, orientations);
% 
% rotatedSeg410 = logical(rotatePharynx(seg410, centroids, orientations));
% rotatedSeg470 = logical(rotatePharynx(seg470, centroids, orientations));

%% Generate Midlines
midlines = calculateMidlines(imTL, seg410);

%% get profile data
i410 = measureIntensityAlongMidlines(im410, midlines);
i470 = measureIntensityAlongMidlines(im470, midlines);

i410_trimmed = trimProfile(i410);
i470_trimmed = trimProfile(i470);

%% Plot Strain Means

% Make a figures directory
mkdir(rawImgDir, 'figures');

plotStrainMeans(strains, i410_trimmed, "Mean 410nm \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_410nm_raw.pdf'));
plotStrainMeans(strains, i470_trimmed, "Mean 470nm \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_470nm_raw.pdf'));

% Calc R, OxD, E
R = i410_trimmed ./ i470_trimmed;
OxD = ja_oxd(R);
E = ja_E(OxD);

plotStrainMeans(strains, E, "Mean E \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_E_raw.pdf'));
plotStrainMeans(strains, OxD, "Mean OxD \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_OxD_raw.pdf'));

%% Save Profile Data
csvwrite(fullfile(rawImgDir, 'i410.csv'),i410_trimmed.');
csvwrite(fullfile(rawImgDir, 'i470.csv'),i470_trimmed.');